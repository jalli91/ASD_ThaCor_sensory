% ABIDE_ASD_preprocess_multisite.m
%
% Multi-site CONN preprocessing script for ABIDE I & II adult ASD subjects.
% Targets subject folders starting with 'ma' (male ASD) or 'fa' (female ASD)
% across all ABIDEI-* and ABIDEII-* institute directories under ImageFiles/.
%
% Folder naming convention applied by rename_subjects.py:
%   ma#### = Male, ASD          fa#### = Female, ASD
%   mc#### = Male, Control      fc#### = Female, Control
%
% Pipeline per subject/session:
%   1. CONN Setup + Preprocessing  (default_mni: realign → coreg → segment
%                                   → normalise to MNI 2mm → 6mm FWHM smooth)
%   2. CONN Denoising              (band-pass 0.01–0.08 Hz, WM/CSF/motion/scrubbing)
%   3. CONN First-level Analysis   (bivariate correlation, HRF weighting)
%   ROIs: Harvard-Oxford atlas (CONN built-in) + default resting-state networks
%
% Reference: https://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html
%            https://fcon_1000.projects.nitrc.org/indi/abide/abide_II.html
%            Parameter PDFs in each ABIDEI-<SITE>/ and ABIDEII-<SITE>_X/ directory
%
% Prerequisites:
%   - SPM12 + CONN toolbox on MATLAB path  (run: addpath spm12; addpath conn)
%   - CONN v21+ recommended for native .nii.gz support; otherwise pre-decompress
%
% Output: one CONN project .mat per subject × session in
%         ImageFiles/ABIDEII-<SITE>_X/batch_output/
%
% Notes on multi-scan subjects:
%   - Multiple anatomical scans (OILH: up to anat_5): always uses anat_1 (T1w)
%   - Multiple rest runs (IP: up to rest_4; OILH: rest_2): all runs included
%     as separate CONN sessions within the same project
%   - Multiple sessions (OILH_2 only, 8 subjects have session_2):
%     processed as independent CONN projects (one per session)

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
% IP — Institut Pasteur and Robert Debré Hospital, Paris, France
% Scanner : Philips Intera Achieva 1.5T
% Source  : func_para_IP.pdf
% TE=45ms is notably longer than other 1.5T sites: deliberate tradeoff for
% enhanced BOLD contrast (T2* of grey matter ~65ms at 1.5T → TE=45ms
% approaches the signal-optimum while BNI uses TE=25ms for SNR).
% SENSE disabled despite Philips hardware (unusual for Philips acquisition).
% IMPORTANT: IP subjects have up to 4 rest runs (rest_1 … rest_4 in session_1).
%            All discovered runs are included as separate CONN sessions.
% ------------------------------------------------------------------
site_params.IP_II = struct( ...
    'scanner',           'Philips Intera Achieva 1.5T', ...
    'field_strength_T',  1.5 , ...%?
    'TR',                2.700, ...
    'TE_s',              0.045, ...   % 45 ms — long for 1.5T; optimises BOLD at cost of SNR
    'flip_angle_deg',    90, ...%?
    'nslices',           32, ...
    'slice_thickness_mm',4.0, ...
    'voxel_mm',          [3.59 3.59 4.0], ...
    'slice_order',       'ascending', ...   % PDF: "Slice scan order = default"; Philips single-shot EPI default is sequential ascending. Differs from BNI which explicitly says "ascend". Verify against NIfTI slice timing header if needed.
    'multiband_factor',  1, ...%?
    'nvolumes',          176, ...     % derived: total scan duration 07:55.2s / TR 2.7s = 175.6 ≈ 176
    'SENSE_factor',      1 ...        % SENSE disabled (confirmed: SENSE = no in PDF)
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
skipped_count    = 0;    % subjects skipped (already complete)

% Whitelist: only these subjects will be processed; all others are skipped.
% Set to {} to disable filtering and process all subjects.
target_subjects = {'fa29613', 'fa29623', 'fa29628'};

% Each rest_* run is preprocessed as a separate CONN project (one project
% per run). The first anatomical image (anat_1) is shared across all runs.

% Enumerate all institute directories (ABIDE I and II)
inst_dirs = [dir(fullfile(imageFilesDir, 'ABIDEI-*')); ...
             dir(fullfile(imageFilesDir, 'ABIDEII-*'))];
inst_dirs = inst_dirs([inst_dirs.isdir]);

% to set starting subject according to previous stop point
start_from = 1;
start_counter = 0;

for inst_idx = 1:length(inst_dirs)

    inst_name = inst_dirs(inst_idx).name;    % e.g. 'ABIDEII-BNI_1'
    inst_path = fullfile(imageFilesDir, inst_name);

    % --- Resolve site code -------------------------------------------
    site_code = '';
    for sm = 1:size(site_map, 1)
        if contains(inst_name, site_map{sm, 1})
            site_code = site_map{sm, 2};
            break;
        end
    end
    if isempty(site_code)
        warning('Unknown site for directory: %s — skipping.', inst_name);
        continue;
    end
    if ~isfield(site_params, site_code)
        fprintf('  Site %s parameters not yet defined — skipping %s.\n', site_code, inst_name);
        continue;
    end

    sp = site_params.(site_code);    % shorthand for this site's parameters

    fprintf('\n=== Site: %s | Scanner: %s | TR=%.3fs | Slices=%s | Order=%s ===\n', ...
            site_code, sp.scanner, sp.TR, ...
            num2str(sp.nslices), sp.slice_order);

    % Warn if multiband — slice timing correction needs attention
    if sp.multiband_factor > 1
        fprintf('    [MB%d] Multiband site: slice_order string is an approximation.\n', ...
                sp.multiband_factor);
        fprintf('    For accurate slice timing correction, extract the timing vector from\n');
        fprintf('    the NIfTI header: hdr=spm_vol(func_file); t=hdr(1).private.timing.toffst\n');
        fprintf('    and set batch.Setup.preprocessing.sliceorder = t;\n');
    end

    % --- Per-site output directory -----------------------------------
    if use_per_site_output
        siteOutDir = fullfile(inst_path, 'batch_output');
    else
        siteOutDir = fullfile(imageFilesDir, 'batch_output');
    end
    if ~exist(siteOutDir, 'dir'), mkdir(siteOutDir); end


    % --- Find all subject folders and order: fa/fc first, then ma/mc ---
    % fa = female ASD, fc = female control, ma = male ASD, mc = male control.
    % Female subjects are processed first so their CONN projects are available
    % for QC and analysis while male subjects are still running.
    all_entries = dir(inst_path);
    all_entries = all_entries([all_entries.isdir]);
    all_names   = {all_entries.name};
    all_mask    = startsWith(all_names, 'fa') | startsWith(all_names, 'fc') | ...
                  startsWith(all_names, 'ma') | startsWith(all_names, 'mc');
    all_subj    = all_names(all_mask);

    % Partition into female-first then male, each group sorted alphabetically
    fa_fc      = sort(all_subj(startsWith(all_subj, 'fa') | startsWith(all_subj, 'fc')));
    ma_mc      = sort(all_subj(startsWith(all_subj, 'ma') | startsWith(all_subj, 'mc')));
    subj_names = [fa_fc, ma_mc];

    fprintf('  Found %d subjects total (%d fa/fc first, then %d ma/mc)\n', ...
            length(subj_names), length(fa_fc), length(ma_mc));

    % ================================================================
    for subj_idx = 1:length(subj_names)

        % start from setting start subject
        start_counter = start_counter + 1;
        if start_counter < start_from
            continue;
        end

        subj_id   = subj_names{subj_idx};

        % Skip subjects not in the whitelist (if whitelist is active)
        if ~isempty(target_subjects) && ~ismember(subj_id, target_subjects)
            continue;
        end
        subj_path = fullfile(inst_path, subj_id);

        fprintf('  [%d/%d] %s\n', subj_idx, length(subj_names), subj_id);

        % --- Discover sessions — use session_1 only (first session if absent) ---
        sess_dirs = dir(fullfile(subj_path, 'session_*'));
        sess_dirs = sess_dirs([sess_dirs.isdir]);

        if isempty(sess_dirs)
            warning('    No session_* folder for %s — skipping.', subj_id);
            failed_subjects{end+1} = subj_id;
            continue;
        end

        % Prefer session_1; fall back to the first available session
        sess_names = sort({sess_dirs.name});
        if ismember('session_1', sess_names)
            sess_dirs = sess_dirs(strcmp({sess_dirs.name}, 'session_1'));
        else
            sess_dirs = sess_dirs(1);
            fprintf('    [WARN] session_1 absent; using %s\n', sess_dirs(1).name);
        end

        for ses_idx = 1:length(sess_dirs)

            ses_name = sess_dirs(ses_idx).name;    % 'session_1' or 'session_2'
            ses_path = fullfile(subj_path, ses_name);

            % ---------------------------------------------------------
            % STRUCTURAL FILE SELECTION
            % ---------------------------------------------------------
            % Always use anat_1 as the T1-weighted structural reference.
            % OILH subjects collected multiple anatomical scans (anat_1 to
            % anat_5/6) — these may include repeat T1, T2, or PDw sequences.
            % Using only anat_1 keeps preprocessing consistent across sites.
            % A fallback to any anat_* is provided if anat_1 is absent.
            anat_path = fullfile(ses_path, 'anat_1');
            if ~exist(anat_path, 'dir')
                anat_candidates = dir(fullfile(ses_path, 'anat_*'));
                anat_candidates = anat_candidates([anat_candidates.isdir]);
                if isempty(anat_candidates)
                    warning('    No anat_* for %s/%s — skipping.', subj_id, ses_name);
                    failed_subjects{end+1} = sprintf('%s_%s', subj_id, ses_name);
                    continue;
                end
                anat_path = fullfile(ses_path, anat_candidates(1).name);
                fprintf('    [WARN] anat_1 absent; using %s\n', anat_candidates(1).name);
            end

            anat_nii = find_nii(anat_path, 'anat');
            if isempty(anat_nii)
                warning('    No NIfTI in %s — skipping.', anat_path);
                failed_subjects{end+1} = sprintf('%s_%s', subj_id, ses_name);
                continue;
            end
            STRUCTURAL_FILE = {anat_nii};

            % ---------------------------------------------------------
            % FUNCTIONAL FILE SELECTION
            % ---------------------------------------------------------
            % Discover all rest_N directories and collect NIfTI paths.
            % CONN treats each rest run as a separate session within the
            % same subject (batch.Setup.functionals{1}{run_idx}{1}).
            %   - BNI/ETH/IU/KUL/NYU/SDSU/TCD/USM: single rest_1 per session
            % Use rest_1 only (first run in folder if rest_1 absent)
            rest_entries = dir(fullfile(ses_path, 'rest_*'));
            rest_entries = rest_entries([rest_entries.isdir]);
            rest_names   = sort({rest_entries.name});   % ascending: rest_1, rest_2, …

            if isempty(rest_names)
                warning('    No rest_* folder for %s/%s — skipping.', subj_id, ses_name);
                failed_subjects{end+1} = sprintf('%s_%s', subj_id, ses_name);
                continue;
            end

            % Discover NIfTI paths for every rest_* run in this session.
            % Each run will become its own CONN project (one project per run).
            rest_run_files = {};
            rest_run_names = {};
            for ri = 1:numel(rest_names)
                func_nii = find_nii(fullfile(ses_path, rest_names{ri}), 'rest');
                if isempty(func_nii)
                    warning('    No NIfTI in %s — skipping this run.', rest_names{ri});
                    continue;
                end
                rest_run_files{end+1} = func_nii;   %#ok<AGROW>
                rest_run_names{end+1} = rest_names{ri}; %#ok<AGROW>
            end

            if isempty(rest_run_files)
                warning('    No valid rest runs for %s/%s — skipping.', subj_id, ses_name);
                failed_subjects{end+1} = sprintf('%s_%s', subj_id, ses_name);
                continue;
            end

            fprintf('    Found %d rest run(s): %s\n', numel(rest_run_files), ...
                    strjoin(rest_run_names, ', '));

            % =============================================================
            % ONE CONN PROJECT PER REST RUN
            % =============================================================
            for ri = 1:numel(rest_run_files)

                run_name        = rest_run_names{ri};         % e.g. 'rest_1'
                FUNCTIONAL_FILES = {rest_run_files{ri}};      % single run

                % ---------------------------------------------------------
                % CONN PROJECT FILE — named per subject/session/run
                % ---------------------------------------------------------
                proj_name = sprintf('%s_%s_%s.mat', subj_id, ses_name, run_name);
                proj_path = fullfile(siteOutDir, proj_name);

                if exist(proj_path, 'file')
                    proj_dir = strrep(proj_path, '.mat', '');
                    if ~isempty(target_subjects) && ismember(subj_id, target_subjects)
                        % Force overwrite for whitelisted subjects
                        fprintf('    [OVERWRITE] Deleting existing output for %s\n', proj_name);
                        delete(proj_path);
                        if exist(proj_dir, 'dir'), rmdir(proj_dir, 's'); end
                    elseif exist(fullfile(proj_dir, 'results', 'firstlevel'), 'dir')
                        fprintf('    [%d done] Already complete, skipping: %s\n', ...
                                skipped_count + completed_count, proj_name);
                        skipped_count = skipped_count + 1;
                        continue;
                    else
                        fprintf('    Incomplete (no first-level), reprocessing: %s\n', proj_name);
                        delete(proj_path);
                        if exist(proj_dir, 'dir'), rmdir(proj_dir, 's'); end
                    end
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
                batch.Setup.structurals = STRUCTURAL_FILE;  % always anat_1

                % Single functional run as one CONN session
                batch.Setup.functionals        = {{{FUNCTIONAL_FILES{1}}}};
                batch.Setup.conditions.names   = {'rest'};
                batch.Setup.conditions.onsets{1}{1}{1}    = 0;
                batch.Setup.conditions.durations{1}{1}{1} = Inf;

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
                    anat_dir_clean = fileparts(STRUCTURAL_FILE{1});

                    % Functional intermediates — clean this run's directory
                    func_prefixes = {'dswau', 'swau', 'wau', 'au', 'u'};
                    func_dir_clean = fileparts(FUNCTIONAL_FILES{1});
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

                    % Anatomical intermediates — only clean after the last run
                    % to avoid removing files still needed by subsequent runs.
                    if ri == numel(rest_run_files)
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
                    end

                catch ME
                    warning('    FAILED: %s\n    Error: %s', proj_name, ME.message);
                    failed_subjects{end+1} = proj_name;
                end

            end   % rest runs loop

        end   % sessions loop
    end       % subjects loop
end           % institutes loop

%% ========================================================
%  SECTION 4 — SUMMARY REPORT
% ========================================================
fprintf('\n\n=== PREPROCESSING RUN COMPLETE ===\n');
fprintf('  Newly completed : %d\n', completed_count);
fprintf('  Already done (skipped) : %d\n', skipped_count);
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
