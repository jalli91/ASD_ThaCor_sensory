% ABIDE_empty_analysis_rescue.m
%
% Rescue preprocessing script for subjects with missing/empty first-level
% analysis results across multiple sites.
%
% Target subjects:
%   ABIDEII-OILH_2 / mc28728    (male Control, OILH site)
%   ABIDEI-Caltech / mc0051488  (male Control, Caltech site)
%
% For each subject, only the FIRST available session and the FIRST
% available anatomical (anat_1) and resting (rest_1) images are used.
%
% Pipeline per subject:
%   1. CONN Setup + Preprocessing  (default_mni: realign → coreg → segment
%                                   → normalise to MNI 2mm → 6mm FWHM smooth)
%   2. CONN Denoising              (band-pass 0.01–0.08 Hz, WM/CSF/motion/scrubbing)
%   3. CONN First-level Analysis   (bivariate correlation, HRF weighting)
%   ROIs: Harvard-Oxford atlas (CONN built-in) + default resting-state networks
%
% Prerequisites:
%   - SPM12 + CONN toolbox on MATLAB path  (run: addpath spm12; addpath conn)
%   - CONN v21+ recommended for native .nii.gz support; otherwise pre-decompress
%
% Output: one CONN project .mat per subject in each site's batch_output/

clear all;

%% ========================================================
%  SECTION 1 — SITE PARAMETER LIBRARY
% ========================================================
% Each struct field corresponds to one ABIDE II site code.
% Key parameters for preprocessing:
%   TR           : repetition time (seconds) — passed to CONN as batch.Setup.RT
%   slice_order  : string passed to batch.Setup.preprocessing.sliceorder
%                  CONN-recognised strings used here:
%                    'ascending'                 — Philips default (bottom→top)
%                    'interleaved (Siemens)'     — Siemens odd-slices-first
%                    'interleaved (AFNI: bottom-up)' — GE interleaved ascending
%   multiband_factor > 1 : simultaneous multi-slice (SMS/MB); slice timing string
%                  is still an approximation — see per-site notes below.

% ------------------------------------------------------------------
% OILH — Olin Neuropsychiatry Research Center, Institute of Living at Hartford Hospital
%         Hartford CT, USA  (ABIDE II code: ABIDEII-ONRC_2; local directory: ABIDEII-OILH_2)
% Scanner : Siemens MAGNETOM Skyra syngo MR D13 3T
% Source  : func_OILH.pdf, anat_OILH.pdf
% *** SIMULTANEOUS MULTI-SLICE (SMS/Multiband), MB factor = 8 ***
% 48 slices with MB8 → 6 simultaneous groups. TR=475ms is the shortest here.
% 947 volumes × 475ms ≈ 7.5 min scan — high temporal resolution.
% Flip angle 60° = Ernst angle for short TR.
% Phase encoding A >> P (only A>>P subjects are preprocessed).
% *** MULTIPLE SESSIONS ***: 8 subjects (fc28711, fc28713, ma28681–28683,
%     ma28687, mc28705, mc28712) have both session_1 and session_2.
%     Each session is processed as an independent CONN project.
% *** MULTIPLE ANATOMICALS ***: up to anat_5 per session; anat_1 used here.
% *** MULTIPLE REST RUNS ***: rest_1 (all subjects) + rest_2 (subset).
%     Both runs are included as separate CONN sessions within one project.
% ------------------------------------------------------------------
site_params.OILH_II = struct( ...
    'scanner',           'Siemens MAGNETOM Skyra syngo MR D13 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                0.475, ...   % very short due to MB8
    'TE_s',              0.030, ...
    'flip_angle_deg',    60, ...
    'nslices',           48, ...
    'slice_thickness_mm',3.0, ...
    'voxel_mm',          [3.0 3.0 3.0], ...
    'slice_order',       'interleaved (Siemens)', ...  % MB8 approx; use timing vector if available
    'multiband_factor',  8, ...       % 48 slices / MB8 = 6 simultaneous groups
    'nvolumes',          947, ...
    'phase_enc_dir',     'A>>P' ...   % only A>>P subjects are preprocessed
);

% ------------------------------------------------------------------
% Caltech_I — California Institute of Technology, Pasadena CA, USA
% Scanner : Siemens MAGNETOM TrioTim syngo MR B17 3T
% Source  : func_Caltech.pdf, anat_Caltech.pdf
% 34 slices (odd count) → Siemens interleaved odd-first.
% No PAT acceleration. 150 volumes × TR 2s = 5 min scan.
% ------------------------------------------------------------------
site_params.Caltech_I = struct( ...
    'scanner',           'Siemens MAGNETOM TrioTim syngo MR B17 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                2.000, ...
    'TE_s',              0.030, ...
    'flip_angle_deg',    75, ...
    'nslices',           34, ...
    'slice_thickness_mm',3.5, ...
    'voxel_mm',          [3.5 3.5 3.5], ...
    'slice_order',       'interleaved (Siemens)', ...   % 34 slices, odd-first
    'multiband_factor',  1, ...
    'nvolumes',          150, ...
    'phase_enc_dir',     'A>>P', ...
    'SENSE_factor',      1 ...
);

%% ========================================================
%  SECTION 2 — DIRECTORY AND PATH CONFIGURATION
% ========================================================

% Root directory containing ABIDEII-*_* institute subdirectories
imageFilesDir = '/Volumes/Milk/CYL/ASD_ThaCor/ABIDE/ImageFiles';

% Output: one CONN project .mat per subject/session, stored inside each
% site's own batch_output/ subfolder for organisation by site.
use_per_site_output = true;

% Map substring in institute directory name → site_params field name.
% ABIDE I entries use 'ABIDEI-<SITE>' prefix to avoid collision with
% ABIDE II sites that share the same short name (NYU, USM).
% ABIDE II entries likewise use 'ABIDEII-' prefix for the ambiguous sites.
% Order matters: more specific patterns should come first.
site_map = { ...
    % --- ABIDE I ---
    'ABIDEI-Caltech',  'Caltech_I'; ...
    'ABIDEI-CMU',      'CMU_I';     ...
    'ABIDEI-MaxMun',   'MaxMun_I';  ...
    'ABIDEI-NYU',      'NYU_I';     ...
    'ABIDEI-Olin',     'Olin_I';    ...
    'ABIDEI-UM_1',     'UM1_I';     ...
    % --- ABIDE II ---
    'ABIDEII-IP',      'IP_II';     ...
    'ABIDEII-IU',      'IU_II';     ...
    'ABIDEII-NYU',     'NYU_II';    ...
    'ABIDEII-OILH',    'OILH_II';   ...
    'ABIDEII-USM',     'USM_II'     ...
};

%% ========================================================
%  SECTION 3 — MAIN PREPROCESSING LOOP
% ========================================================

failed_subjects  = {};   % accumulate failures for end-of-run report
completed_count  = 0;    % subjects fully processed this run

% --- Rescue target list: {inst_name, subj_id, site_params_key} -------
% Each row identifies one subject to (re)process.
rescue_list = { ...
    'ABIDEII-OILH_2',  'mc28728',   'OILH_II';   ...
    'ABIDEI-Caltech',  'mc0051488', 'Caltech_I'   ...
};

fprintf('  Rescue list: %d subjects\n', size(rescue_list, 1));

for subj_idx = 1:size(rescue_list, 1)

    inst_name = rescue_list{subj_idx, 1};
    subj_id   = rescue_list{subj_idx, 2};
    sp_key    = rescue_list{subj_idx, 3};
    sp        = site_params.(sp_key);

    inst_path = fullfile(imageFilesDir, inst_name);

    if use_per_site_output
        siteOutDir = fullfile(inst_path, 'batch_output');
    else
        siteOutDir = fullfile(imageFilesDir, 'batch_output');
    end
    if ~exist(siteOutDir, 'dir'), mkdir(siteOutDir); end

    fprintf('\n=== Rescue: %s | Scanner: %s | TR=%.3fs | Slices=%d | Order=%s ===\n', ...
            inst_name, sp.scanner, sp.TR, sp.nslices, sp.slice_order);
    if sp.multiband_factor > 1
        fprintf('    [MB%d] Multiband site: slice_order string is an approximation.\n', ...
                sp.multiband_factor);
    end

    subj_path = fullfile(inst_path, subj_id);

    fprintf('  [%d/%d] %s\n', subj_idx, size(rescue_list, 1), subj_id);

    % ---------------------------------------------------------
    % STRUCTURAL FILE SELECTION
    % ---------------------------------------------------------
    % These subjects store all anatomicals in one session and all
    % rest runs in a separate session.  Search every session_* for
    % the first anat_* directory, independent of which session it is in.
    sess_dirs = dir(fullfile(subj_path, 'session_*'));
    sess_dirs = sess_dirs([sess_dirs.isdir]);
    sess_names_sorted = sort({sess_dirs.name});

    anat_nii = '';
    anat_ses_used = '';
    for sn = sess_names_sorted
        anat_path = fullfile(subj_path, sn{1}, 'anat_1');
        if ~exist(anat_path, 'dir')
            anat_cands = dir(fullfile(subj_path, sn{1}, 'anat_*'));
            anat_cands = anat_cands([anat_cands.isdir]);
            if isempty(anat_cands), continue; end
            anat_path = fullfile(subj_path, sn{1}, anat_cands(1).name);
            fprintf('    [WARN] anat_1 absent in %s; using %s\n', sn{1}, anat_cands(1).name);
        end
        anat_nii = find_nii(anat_path, 'anat');
        if ~isempty(anat_nii)
            anat_ses_used = sn{1};
            break;
        end
    end

    if isempty(anat_nii)
        warning('    No anat NIfTI found for %s — skipping.', subj_id);
        failed_subjects{end+1} = subj_id;
        continue;
    end
    STRUCTURAL_FILE = {anat_nii};
    fprintf('    Anat  : %s\n', anat_nii);

    % ---------------------------------------------------------
    % FUNCTIONAL FILE SELECTION
    % ---------------------------------------------------------
    % Search every session_* for the first rest_* directory,
    % independent of which session contains the anatomical.
    func_nii = '';
    rest_ses_used = '';
    for sn = sess_names_sorted
        rest_entries = dir(fullfile(subj_path, sn{1}, 'rest_*'));
        rest_entries = rest_entries([rest_entries.isdir]);
        rest_names   = sort({rest_entries.name});
        if isempty(rest_names), continue; end

        if ismember('rest_1', rest_names)
            selected_rest = 'rest_1';
        else
            selected_rest = rest_names{1};
            fprintf('    [WARN] rest_1 absent in %s; using %s\n', sn{1}, selected_rest);
        end

        func_nii = find_nii(fullfile(subj_path, sn{1}, selected_rest), 'rest');
        if ~isempty(func_nii)
            rest_ses_used = sn{1};
            break;
        end
    end

    if isempty(func_nii)
        warning('    No rest NIfTI found for %s — skipping.', subj_id);
        failed_subjects{end+1} = subj_id;
        continue;
    end
    FUNCTIONAL_FILES = {func_nii};
    nruns = 1;
    fprintf('    Rest  : %s\n', func_nii);

    % ---------------------------------------------------------
    % CONN PROJECT FILE
    % ---------------------------------------------------------
    % Name the project after the subject + the session that contains the
    % resting file (e.g. mc28712_session_2.mat), so multiple sessions
    % produce distinct CONN projects in batch_output/.
    proj_name = sprintf('%s_%s.mat', subj_id, rest_ses_used);
    proj_path = fullfile(siteOutDir, proj_name);

    if exist(proj_path, 'file')
        proj_dir = strrep(proj_path, '.mat', '');
        fprintf('    Overwriting existing project: %s\n', proj_name);
        delete(proj_path);
        if exist(proj_dir, 'dir'), rmdir(proj_dir, 's'); end
    end

    % ---------------------------------------------------------
    % BUILD CONN BATCH STRUCTURE
    % ---------------------------------------------------------
    clear batch;
    batch.filename = proj_path;

    %% --- Setup & Preprocessing ---
    batch.Setup.isnew      = 1;
    batch.Setup.nsubjects  = 1;
    batch.Setup.RT         = sp.TR;         % TR (seconds)
    batch.Setup.structurals = STRUCTURAL_FILE;

    % Assign functional runs as CONN sessions (one file per session)
    batch.Setup.functionals = repmat({{}}, [1, 1]);
    for ri = 1:nruns
        batch.Setup.functionals{1}{ri}{1} = FUNCTIONAL_FILES{ri};
    end

    % Resting-state condition: onset=0, duration=Inf for every run
    batch.Setup.conditions.names = {'rest'};
    for ri = 1:nruns
        batch.Setup.conditions.onsets{1}{1}{ri}    = 0;
        batch.Setup.conditions.durations{1}{1}{ri} = Inf;
    end

    % Preprocessing pipeline: standard MNI normalisation
    %   Realignment (motion correction) → coregistration (func→anat) →
    %   segmentation (GM/WM/CSF) → normalisation (MNI 2mm isotropic) →
    %   spatial smoothing (6mm FWHM Gaussian)
    % Site-specific TR and slice order drive slice timing correction.
    batch.Setup.preprocessing.steps      = 'default_mni';
    batch.Setup.preprocessing.fwhm       = 6;             % smoothing kernel mm
    batch.Setup.preprocessing.sliceorder = sp.slice_order;

    %% --- ROI / Atlas ---
    % Harvard-Oxford cortical atlas (CONN built-in atlas.nii) +
    % Default resting-state networks (networks.nii)
    % These match the settings in ABIDE_preprocess_runner.m
    batch.Setup.rois.names   = {'HO_atlas', 'Default_network'};
    batch.Setup.rois.files{1} = fullfile(fileparts(which('conn')), 'rois', 'atlas.nii');
    batch.Setup.rois.files{2} = fullfile(fileparts(which('conn')), 'rois', 'networks.nii');

    batch.Setup.done      = 1;
    batch.Setup.overwrite = 'Yes';

    %% --- Denoising ---
    % Band-pass filter: 0.01–0.08 Hz (standard BOLD resting-state band).
    % For MB sites (IU TR=813ms, OILH TR=475ms), the Nyquist frequency
    % is 0.61/1.05 Hz — consider expanding the high-pass cutoff to 0.1Hz
    % or using ICA-FIX/aCompCor for physiological noise removal.
    % Confound regressors:
    %   White Matter (6 PCs)  — non-neural BOLD fluctuations
    %   CSF          (6 PCs)  — cardiac/respiratory noise proxy
    %   Realignment  (6 PCs)  — 3 translations + 3 rotations
    %   Scrubbing    (Inf)    — spike regressors for high-motion volumes
    %                           (FD > 0.5mm or DVARS > 3 by default)
    %   Effect of rest (Inf)  — mean-centre each session
    batch.Denoising.filter               = [0.01, 0.08];
    batch.Denoising.confounds.names      = {'White Matter', 'CSF', ...
                                            'realignment', 'scrubbing', ...
                                            'Effect of rest'};
    batch.Denoising.confounds.dimensions = {6 6 6 Inf Inf};
    batch.Denoising.done                 = 1;
    batch.Denoising.overwrite            = 'Yes';

    %% --- First-level Analysis ---
    batch.Analysis.analysis_number = 1;
    batch.Analysis.measure         = 1;   % bivariate correlation
    batch.Analysis.weight          = 2;   % HRF weighting
    batch.Analysis.sources         = {};  % all defined ROIs as seeds
    batch.Analysis.done            = 1;
    batch.Analysis.overwrite       = 'Yes';

    %% --- Run ---
    try
        conn_batch(batch);
        completed_count = completed_count + 1;
        fprintf('    Done [%d complete so far]: %s\n', completed_count, proj_name);

        % --- Delete SPM/CONN intermediate images ---
        % Prefixes ordered longest-first so more specific patterns
        % are reported before the shorter overlapping ones.
        func_dir_clean = fileparts(FUNCTIONAL_FILES{1});
        anat_dir_clean = fileparts(STRUCTURAL_FILE{1});

        % Functional intermediates (realign/coreg/norm/smooth prefixes)
        func_prefixes = {'dswau', 'swau', 'wau', 'au', 'u'};
        for pi = 1:numel(func_prefixes)
            pfx = func_prefixes{pi};
            for ext = {'*.nii', '*.nii.gz'}
                hits = dir(fullfile(func_dir_clean, [pfx ext{1}]));
                for hi = 1:numel(hits)
                    fpath = fullfile(hits(hi).folder, hits(hi).name);
                    delete(fpath);
                    fprintf('      Deleted functional intermediate: %s\n', hits(hi).name);
                end
            end
        end

        % Anatomical intermediates (segmentation/normalisation outputs)
        anat_prefixes = {'ewc', 'iy_c', 'wc'};
        for pi = 1:numel(anat_prefixes)
            pfx = anat_prefixes{pi};
            for ext = {'*.nii', '*.nii.gz'}
                hits = dir(fullfile(anat_dir_clean, [pfx ext{1}]));
                for hi = 1:numel(hits)
                    fpath = fullfile(hits(hi).folder, hits(hi).name);
                    delete(fpath);
                    fprintf('      Deleted anatomical intermediate: %s\n', hits(hi).name);
                end
            end
        end

    catch ME
        warning('    FAILED: %s\n    Error: %s', proj_name, ME.message);
        failed_subjects{end+1} = proj_name;
    end

end       % subjects loop

%% ========================================================
%  SECTION 4 — SUMMARY REPORT
% ========================================================
fprintf('\n\n=== PREPROCESSING RUN COMPLETE ===\n');
fprintf('  Newly completed : %d\n', completed_count);
fprintf('  Failed          : %d\n', length(failed_subjects));
if ~isempty(failed_subjects)
    fprintf('  Failed subjects:\n');
    for fi = 1:length(failed_subjects)
        fprintf('    - %s\n', failed_subjects{fi});
    end
    fprintf('Review warnings above for details.\n');
end

%% ========================================================
%  HELPER FUNCTION
% ========================================================
function nii_path = find_nii(search_dir, prefix)
% find_nii  Locate the primary NIfTI file for a given modality prefix.
%
%   Searches search_dir for <prefix>.nii (uncompressed) then <prefix>.nii.gz.
%   Falls back to the first *.nii / *.nii.gz file if an exact prefix match
%   is not found (handles sites that name files differently).
%
%   If only a .nii.gz is found, it is decompressed to .nii in-place so that
%   SPM/CONN can read it directly.  If the .nii already exists (from a prior
%   run), gunzip is skipped — the existing file is used as-is.
%
%   Returns: full file path string, or '' if nothing found.

    % Try exact name matches first
    candidates = [
        dir(fullfile(search_dir, [prefix '.nii']));
        dir(fullfile(search_dir, [prefix '.nii.gz']))
    ];

    % Fallback: any NIfTI file in the directory
    if isempty(candidates)
        candidates = [
            dir(fullfile(search_dir, '*.nii'));
            dir(fullfile(search_dir, '*.nii.gz'))
        ];
        % Exclude 4D companion JSON / header files if present
        candidates = candidates(~contains({candidates.name}, 'json'));
    end

    % Exclude macOS resource-fork files (._*) — they are AppleDouble metadata,
    % not real NIfTI/GZIP files, and will cause gunzip to fail.
    if ~isempty(candidates)
        candidates = candidates(~startsWith({candidates.name}, '._'));
    end

    if isempty(candidates)
        nii_path = '';
        return;
    end

    % Prefer uncompressed .nii over .nii.gz (lower I/O overhead for SPM)
    names = {candidates.name};
    uncomp = names(endsWith(names, '.nii') & ~endsWith(names, '.nii.gz'));
    if ~isempty(uncomp)
        nii_path = fullfile(search_dir, uncomp{1});
    else
        % Only a .nii.gz is available — decompress if .nii not already present
        gz_path  = fullfile(search_dir, names{1});
        nii_path = strrep(gz_path, '.gz', '');
        if exist(nii_path, 'file')
            fprintf('      [skip gunzip] Already unzipped: %s\n', ...
                    nii_path);
        else
            fprintf('      Decompressing: %s\n', names{1});
            gunzip(gz_path, search_dir);
        end
    end
end
