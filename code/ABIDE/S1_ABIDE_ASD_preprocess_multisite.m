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
% BNI — Barrow Neurological Institute, Phoenix AZ, USA
% Scanner : Philips Intera Achieva 1.5T
% Source  : func_para_BNI.pdf
% TR=3s is the longest in this dataset, giving lower temporal resolution.
% SENSE×2 (AP) reduces acquisition time; 50 axial slices, 4mm thick, no gap.
% Ascending order is the Philips factory default for single-shot EPI.
% ------------------------------------------------------------------
% site_params.BNI_II = struct( ...
%     'scanner',           'Philips Intera Achieva 1.5T', ...
%     'field_strength_T',  1.5, ...
%     'TR',                3.000, ...   % seconds
%     'TE_s',              0.025, ...   % 25 ms
%     'flip_angle_deg',    80, ...
%     'nslices',           50, ...
%     'slice_thickness_mm',4.0, ...
%     'voxel_mm',          [3.75 3.75 4.0], ...
%     'slice_order',       'ascending', ...
%     'multiband_factor',  1, ...
%     'nvolumes',          120, ...
%     'SENSE_factor',      2 ...
% );

% ------------------------------------------------------------------
% ETH — ETH Zurich, Switzerland
% Source  : func_para_ETH.pdf (confirmed); anat_para_ETH.pdf
% Scanner note: anat_para_ETH.pdf lists Siemens MAGNETOM Allegra syngo 3T,
%   but func_para_ETH.pdf uses Philips FFE/SENSE/SPIR parameter format —
%   ETH acquired structural on Siemens Allegra and functional on a Philips
%   scanner (two-scanner site). Functional parameters below are from the
%   Philips functional protocol sheet.
% Slice scan order = "descend" (top→bottom, Philips 'descending').
%   CONN string: 'descending'
% SENSE×2.5 (AP direction); 40 transverse slices; 0.3mm inter-slice gap.
% 210 volumes × TR 2s = 7:06 total scan duration.
% TE=25ms at 3T is short (T2* cortex ~35ms at 3T): trades some BOLD
% sensitivity for reduced susceptibility-induced geometric distortion.
% ------------------------------------------------------------------
% site_params.ETH_II = struct( ...
%     'scanner',           'Philips (functional) / Siemens Allegra 3T (structural)', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.000, ...   % confirmed: TR(ms) = 2000
%     'TE_s',              0.025, ...   % confirmed: TE(ms) = 25
%     'flip_angle_deg',    90, ...      % confirmed
%     'nslices',           40, ...      % confirmed: slices = 40
%     'slice_thickness_mm',3.0, ...     % confirmed: slice thickness = 3mm
%     'slice_gap_mm',      0.3, ...     % confirmed: gap(mm) = 0.299999952
%     'voxel_mm',          [3.0 3.0 3.0], ...  % REC voxel MPS: 3.00/3.00/3.00
%     'slice_order',       'descending', ...    % confirmed: Slice scan order = "descend"
%     'multiband_factor',  1, ...
%     'nvolumes',          210, ...     % confirmed: dyn scans = 210
%     'SENSE_factor',      2.5 ...      % confirmed: P reduction (AP) = 2.5
% );

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

% ------------------------------------------------------------------
% IU — Indiana University, Bloomington IN, USA
% Scanner : Siemens MAGNETOM TrioTim syngo MR B17 3T
% Source  : rest_para_IU.pdf
% *** SIMULTANEOUS MULTI-SLICE (SMS/Multiband), MB factor = 3 ***
% 42 slices acquired in simultaneous groups of 3 (14 effective timepoints/TR).
% Very short TR (813 ms) enables 1200 volumes in ~16 min — highest temporal
% density in this dataset. Flip angle 60° = Ernst angle for TR~1s BOLD.
% For accurate slice timing correction with MB acquisitions CONN/SPM need the
% actual slice onset vector (stored in the NIfTI slice_code/timing field if
% converted with dcm2niix). A fallback interleaved string is set here, but
% consider extracting: hdr = spm_vol(func); t = hdr(1).private.timing.toffst
% and passing it as a numeric vector to batch.Setup.preprocessing.sliceorder.
% Phase encoding: A >> P
% ------------------------------------------------------------------
site_params.IU_II = struct( ...
    'scanner',           'Siemens MAGNETOM TrioTim syngo MR B17 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                0.813, ...
    'TE_s',              0.028, ...
    'flip_angle_deg',    60, ...      % Ernst angle approx for short TR
    'nslices',           42, ...
    'slice_thickness_mm',3.4, ...
    'voxel_mm',          [3.4 3.4 3.4], ...
    'slice_order',       'interleaved (Siemens)', ...  % approx; provide timing vector for accuracy
    'multiband_factor',  3, ...       % 42 slices / MB3 = 14 simultaneous groups
    'nvolumes',          1200, ...
    'phase_enc_dir',     'A>>P', ...
    'SENSE_factor',      1 ...
);

% ------------------------------------------------------------------
% KUL — Katholieke Universiteit Leuven, Belgium
% Scanner : Siemens MAGNETOM Allegra syngo MR 2004A 3T
% Source  : func_para_KUL.pdf
% TE=15ms is unusually short (well below T2* ~35ms at 3T for cortex):
% reduces BOLD sensitivity but also minimises susceptibility-induced signal
% drop-out near orbitofrontal cortex / inferior temporal sulcus.
% Phase encoding R >> L is non-standard (most sites use A >> P or L >> R);
% matters for B0 fieldmap-based distortion correction if applied later.
% Only male ASD (ma) subjects from this site in the current dataset.
% 180 volumes × TR 2s = 6 min scan.
% ------------------------------------------------------------------
% site_params.KUL_II = struct( ...
%     'scanner',           'Siemens MAGNETOM Allegra syngo MR 2004A 3T', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.000, ...
%     'TE_s',              0.015, ...   % 15 ms — very short TE, minimal BOLD contrast
%     'flip_angle_deg',    90, ...
%     'nslices',           33, ...
%     'slice_thickness_mm',4.0, ...
%     'voxel_mm',          [3.0 3.0 4.0], ...
%     'slice_order',       'interleaved (Siemens)', ...
%     'multiband_factor',  1, ...
%     'nvolumes',          180, ...
%     'phase_enc_dir',     'R>>L' ...   % non-standard PE direction
% );

% ------------------------------------------------------------------
% NYU — New York University Langone Medical Center, NY, USA
% Scanner : Siemens MAGNETOM Allegra syngo MR 2004A 3T
% Source  : func_NYU1.pdf, anat_NYU1.pdf
% Only 4 subjects in this dataset (ma29203, ma29205, fa29212, fc29242).
% 33 slices (odd) → Siemens interleaved odd-first.
% Phase encoding R >> L (non-standard; same scanner/protocol as ABIDE I NYU_I).
% TE=15ms is unusually short (well below T2* ~35ms at 3T): minimises
% susceptibility-induced signal drop-out near OFC.
% 180 volumes × TR 2s = 6 min scan.
% ------------------------------------------------------------------
site_params.NYU_II = struct( ...
    'scanner',           'Siemens MAGNETOM Allegra syngo MR 2004A 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                2.000, ...
    'TE_s',              0.015, ...
    'flip_angle_deg',    90, ...
    'nslices',           33, ...
    'slice_thickness_mm',4.0, ...
    'voxel_mm',          [3.0 3.0 4.0], ...
    'slice_order',       'interleaved (Siemens)', ...  % odd-first for 33 slices
    'multiband_factor',  1, ...
    'nvolumes',          180, ...
    'phase_enc_dir',     'R>>L' ...   % confirmed from PDF; non-standard PE direction
);

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
% SDSU — San Diego State University, San Diego CA, USA (SDSU_1)
% Scanner : GE (identical protocol to NYU; same PDF content)
% Source  : rest_para_SDSU.pdf
% Only 1 ASD subject in this dataset: fa28871 (female ASD).
% ------------------------------------------------------------------
% site_params.SDSU_II = struct( ...
%     'scanner',           'GE 3T (model unspecified)', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.000, ...
%     'TE_s',              0.030, ...
%     'flip_angle_deg',    90, ...
%     'nslices',           NaN, ...
%     'slice_thickness_mm',3.4, ...
%     'voxel_mm',          [3.4375 3.4375 3.4], ...
%     'slice_order',       'interleaved (AFNI: bottom-up)', ...  % GE bottom-up
%     'multiband_factor',  1, ...
%     'nvolumes',          NaN, ...
%     'dummy_volumes',     5 ...
% );

% ------------------------------------------------------------------
% TCD — Trinity Centre for Health Sciences, Dublin, Ireland
% Scanner : Siemens MAGNETOM TrioTim syngo MR B17 3T
%           (same scanner model as IU)
% Source  : func_para_TCD.pdf
% GRAPPA×2 in PE direction (A >> P) reduces g-factor noise with full-FOV.
% 10% slice gap (dist_factor=10%) between 40 slices — effective slice spacing
% = 3.0 × 1.10 = 3.3 mm. Enter effective thickness for SPM slice timing.
% 240 volumes × TR 2s = 8 min — longest standard-TR acquisition in dataset.
% ------------------------------------------------------------------
% site_params.TCD_II = struct( ...
%     'scanner',           'Siemens MAGNETOM TrioTim syngo MR B17 3T', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.000, ...
%     'TE_s',              0.028, ...
%     'flip_angle_deg',    90, ...
%     'nslices',           40, ...
%     'slice_thickness_mm',3.0, ...
%     'slice_gap_pct',     10, ...      % 10% gap → effective spacing = 3.0×1.10 = 3.3 mm
%     'voxel_mm',          [3.4 3.4 3.0], ...
%     'slice_order',       'interleaved (Siemens)', ...
%     'multiband_factor',  1, ...
%     'nvolumes',          240, ...
%     'GRAPPA_factor',     2, ...
%     'phase_enc_dir',     'A>>P' ...
% );

% ------------------------------------------------------------------
% USM — University of Utah School of Medicine, Salt Lake City UT, USA
% Scanner : Philips (confirmed from rest_para_USM.pdf: FFE/SENSE/SPIR format)
% Source  : rest_para_USM.pdf
% Highest in-plane resolution (2.5mm) of any site here.
% SENSE×2 (AP direction), ascending slice order, SPIR fat suppression.
% 162 dynamic scans × TR 2.5s ≈ 6.75 min.
% ------------------------------------------------------------------
site_params.USM_II = struct( ...
    'scanner',           'Siemens MAGNETOM TrioTim syngo MR B17 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                2.000, ...       % TR 2000 ms
    'TE_s',              0.028, ...       % TE 28 ms
    'flip_angle_deg',    90, ...          % ✓ 吻合
    'nslices',           40, ...          % Slices 40
    'slice_thickness_mm',3.0, ...         % Slice thickness 3.0 mm
    'voxel_mm',          [3.4 3.4 3.0], ...  % Voxel size: 3.4×3.4×3.0 (FOV 220/base 64 = 3.4375 ≈ 3.4)
    'slice_order',       'interleaved (Siemens)', ...  % Interleaved, 40=偶數→偶先
    'multiband_factor',  1, ...           % 無 MB 欄位
    'nvolumes',          240, ...         % Measurements 240; 240×2s=480s=8:00 ≈ TA 8:06
    'slice_gap_mm',      0.3, ...         % Dist. factor 10% × 3.0mm thickness = 0.3mm gap
    'phase_enc_dir',     'A>>P', ...
    'SENSE_factor',      2 ...            % PAT mode GRAPPA, Accel. factor PE = 2
);

%% ========================================================
%  ABIDE I SITE PARAMETERS
% ========================================================
% Scanner parameters extracted from func_<SITE>.pdf and anat_<SITE>.pdf
% stored in each ABIDEI-<SITE>/ directory.
% Slice order for Siemens interleaved: odd slices first (1,3,5…) when
% nslices is odd; even slices first (2,4,6…) when nslices is even.
% Slice order for Philips EPI: sequential ascending (default for Philips).

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

% ------------------------------------------------------------------
% CMU_I — Carnegie Mellon University, Pittsburgh PA, USA
% Scanner : Siemens MAGNETOM Verio syngo MR B17 3T
% Source  : func_CMU.pdf, anat_CMU.pdf
% 28 slices (even count) → Siemens interleaved even-first.
% Dist. factor 50% — 1.5mm gap between 3.0mm slices (effective spacing = 4.5mm).
% No PAT acceleration. 240 volumes × TR 2s = 8 min scan.
% ------------------------------------------------------------------
site_params.CMU_I = struct( ...
    'scanner',           'Siemens MAGNETOM Verio syngo MR B17 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                2.000, ...
    'TE_s',              0.030, ...
    'flip_angle_deg',    73, ...
    'nslices',           28, ...
    'slice_thickness_mm',3.0, ...
    'voxel_mm',          [3.0 3.0 3.0], ...
    'slice_order',       'interleaved (Siemens)', ...   % 28 slices, even-first
    'multiband_factor',  1, ...
    'nvolumes',          240, ...
    'phase_enc_dir',     'A>>P', ...
    'SENSE_factor',      1 ...
);

% ------------------------------------------------------------------
% Leuven_I — Katholieke Universiteit Leuven, Belgium (ABIDE I)
% Scanner : Philips Intera 3T
% Source  : Func_Leuven_1.pdf, anat_Leuven_1.pdf
% TR=1667ms (act.) — short TR for Philips EPI.
% SENSE×2 (AP); 32 slices, 4.0mm thick, 3.59mm in-plane.
% Philips EPI default slice order: sequential ascending.
% 250 dynamic scans × 1.667s ≈ 6.94 min.
% Note: Leuven_1 only has male subjects in this dataset.
% ------------------------------------------------------------------
% site_params.Leuven_I = struct( ...
%     'scanner',           'Philips Intera 3T', ...
%     'field_strength_T',  3.0, ...
%     'TR',                1.667, ...
%     'TE_s',              0.033, ...
%     'flip_angle_deg',    90, ...
%     'nslices',           32, ...
%     'slice_thickness_mm',4.0, ...         % corrected from PDF (was 2.72)
%     'voxel_mm',          [3.59 3.59 4.0], ...  % corrected from PDF (was [2.75 2.75 2.72])
%     'slice_order',       'ascending', ...    % Philips default
%     'multiband_factor',  1, ...
%     'nvolumes',          250, ...
%     'phase_enc_dir',     'A>>P', ...
%     'SENSE_factor',      2 ...            % corrected from PDF (was 3)
% );

% ------------------------------------------------------------------
% MaxMun_I — Max Planck Institute, Munich, Germany (ABIDE I)
% Scanner : Siemens MAGNETOM Verio syngo MR B17 3T
% Source  : func_MaxMun.pdf, anat_MaxMun.pdf
% MaxMun used three rest protocols (Rest_a/b/c); adult subjects used
% Rest_a: TR=3000ms, 28 slices, 4.0mm thick, 3×3×4mm voxels, 120 vols.
% Rest_b (35 slices) and Rest_c (40 slices, 3mm, 200 vols) apply to
% different subject groups — see func_MaxMun.pdf for details.
% 28 slices (even) → Siemens interleaved even-first.
% Longest TR in dataset (3s): lower temporal SNR but good brain coverage.
% ------------------------------------------------------------------
site_params.MaxMun_I = struct( ...
    'scanner',           'Siemens MAGNETOM Verio syngo MR B17 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                3.000, ...   % Rest_a (adult protocol)
    'TE_s',              0.030, ...
    'flip_angle_deg',    80, ...
    'nslices',           28, ...      % Rest_a: 28 slices
    'slice_thickness_mm',4.0, ...     % Rest_a: 4mm thick
    'voxel_mm',          [3.0 3.0 4.0], ...
    'slice_order',       'interleaved (Siemens)', ...   % 28 slices, even-first
    'multiband_factor',  1, ...
    'nvolumes',          120, ...     % Rest_a: 120 volumes
    'phase_enc_dir',     'A>>P', ...
    'SENSE_factor',      1 ...
);

% ------------------------------------------------------------------
% NYU_I — New York University Langone Medical Center, NY, USA (ABIDE I)
% Scanner : Siemens MAGNETOM Allegra syngo MR 2004A 3T
% Source  : func_NYU.pdf, anat_NYU.pdf
% NOTE: Phase encoding R >> L — non-standard direction.
%   This differs from ABIDE II NYU_1 (different scanner & parameters).
% TE=15ms is short (well below cortical T2*~35ms at 3T); reduces BOLD
% sensitivity but minimises distortion (same rationale as KUL in ABIDE II).
% 33 slices (odd) → Siemens interleaved odd-first.
% 180 volumes × TR 2s = 6 min scan.
% ------------------------------------------------------------------
site_params.NYU_I = struct( ...
    'scanner',           'Siemens MAGNETOM Allegra syngo MR 2004A 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                2.000, ...
    'TE_s',              0.015, ...   % 15ms — short TE, minimal distortion
    'flip_angle_deg',    90, ...
    'nslices',           33, ...
    'slice_thickness_mm',4.0, ...
    'voxel_mm',          [3.0 3.0 4.0], ...
    'slice_order',       'interleaved (Siemens)', ...   % 33 slices, odd-first
    'multiband_factor',  1, ...
    'nvolumes',          180, ...
    'phase_enc_dir',     'R>>L', ...  % non-standard — verify distortion correction direction
    'SENSE_factor',      1 ...
);

% ------------------------------------------------------------------
% Olin_I — Olin Neuropsychiatry Research Center, Hartford CT, USA (ABIDE I)
% Scanner : Siemens MAGNETOM Allegra syngo MR 2004A 3T
% Source  : func_Olin.pdf, anat_Olin.pdf
% Short TR (1500ms) provides higher temporal resolution.
% Flip angle 60° ≈ Ernst angle for TR~1.5s BOLD.
% 29 slices (odd) → Siemens interleaved odd-first.
% 210 volumes × 1.5s = 5.25 min scan.
% ------------------------------------------------------------------
site_params.Olin_I = struct( ...
    'scanner',           'Siemens MAGNETOM Allegra syngo MR 2004A 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                1.500, ...
    'TE_s',              0.027, ...
    'flip_angle_deg',    60, ...      % Ernst angle for short TR
    'nslices',           29, ...
    'slice_thickness_mm',4.0, ...
    'voxel_mm',          [3.4 3.4 4.0], ...
    'slice_order',       'interleaved (Siemens)', ...   % 29 slices, odd-first
    'multiband_factor',  1, ...
    'nvolumes',          210, ...
    'phase_enc_dir',     'A>>P', ...
    'SENSE_factor',      1 ...
);

% ------------------------------------------------------------------
% Pitt_I — University of Pittsburgh, Pittsburgh PA, USA (ABIDE I)
% Scanner : Siemens MAGNETOM Allegra syngo MR 2004A 3T
% Source  : func_Pitt.pdf, anat_Pitt.pdf
% Short TR (1500ms), same scanner family as Olin.
% 29 slices (odd) → Siemens interleaved odd-first.
% 200 volumes × 1.5s = 5 min scan.
% ------------------------------------------------------------------
% site_params.Pitt_I = struct( ...
%     'scanner',           'Siemens MAGNETOM Allegra syngo MR 2004A 3T', ...
%     'field_strength_T',  3.0, ...
%     'TR',                1.500, ...
%     'TE_s',              0.025, ...
%     'flip_angle_deg',    70, ...
%     'nslices',           29, ...
%     'slice_thickness_mm',4.0, ...
%     'voxel_mm',          [3.1 3.1 4.0], ...
%     'slice_order',       'interleaved (Siemens)', ...   % 29 slices, odd-first
%     'multiband_factor',  1, ...
%     'nvolumes',          200, ...
%     'phase_enc_dir',     'A>>P', ...
%     'SENSE_factor',      1 ...
% );

% ------------------------------------------------------------------
% SBL_I — Social Brain Lab, Netherlands Institute for Neuroscience,
%          Amsterdam, Netherlands (ABIDE I)
% Scanner : Philips Intera 3T
% Source  : func_SBL.pdf, anat_SBL.pdf
% SENSE×3 (AP); 38 slices, 2.72mm thick, 2.75mm in-plane.
% 200 dynamic scans × TR 2.2s ≈ 7.33 min.
% Slice scan order: descending (confirmed from PDF: Slice scan order = descend).
% ------------------------------------------------------------------
% site_params.SBL_I = struct( ...
%     'scanner',           'Philips Intera 3T', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.200, ...
%     'TE_s',              0.030, ...
%     'flip_angle_deg',    80, ...
%     'nslices',           38, ...
%     'slice_thickness_mm',2.72, ...
%     'voxel_mm',          [2.75 2.75 2.72], ...
%     'slice_order',       'descending', ...   % confirmed from PDF (was 'ascending')
%     'multiband_factor',  1, ...
%     'nvolumes',          200, ...
%     'phase_enc_dir',     'A>>P', ...
%     'SENSE_factor',      3 ...
% );

% ------------------------------------------------------------------
% Trinity_I — Trinity Centre for Health Sciences, Dublin, Ireland (ABIDE I)
% Scanner : Philips Intera 3T
% Source  : Func_Trinity.pdf, anat_Trinity.pdf
% SENSE×2 (AP), 38 slices, 3.5mm thick, 3.0mm in-plane.
% 150 dynamic scans × TR 2s = 5 min scan.
% Philips EPI default slice order: sequential ascending.
% Note: Trinity used Philips Intera 3T; UM_2 used a GE scanner with reverse
%       spiral sequence — protocols differ (see UM2_I below).
% ------------------------------------------------------------------
% site_params.Trinity_I = struct( ...
%     'scanner',           'Philips Intera 3T', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.000, ...
%     'TE_s',              0.028, ...
%     'flip_angle_deg',    90, ...
%     'nslices',           38, ...
%     'slice_thickness_mm',3.5, ...
%     'voxel_mm',          [3.0 3.0 3.5], ...
%     'slice_order',       'ascending', ...    % Philips default
%     'multiband_factor',  1, ...
%     'nvolumes',          150, ...
%     'phase_enc_dir',     'A>>P', ...
%     'SENSE_factor',      2 ...
% );

% ------------------------------------------------------------------
% UM1_I — University of Michigan, Ann Arbor MI, USA — Scanner 1 (ABIDE I)
% Scanner : Siemens MAGNETOM Allegra syngo MR A30 3T
% Source  : func_UM_1.pdf, anat_UM_1.pdf
% Protocol is ep2d_bold — standard Cartesian EPI (not spiral).
% Parameters are nearly identical to Pitt_I (same Allegra A30 scanner family,
%   same TR/TE/FA/slices/thickness/voxels/PAT).
% 29 slices (odd) → Siemens interleaved odd-first.
% 200 volumes × TR 1.5s = 5:00 scan (Scan Time confirmed: 5:06 in PDF).
% No PAT acceleration.
% ------------------------------------------------------------------
site_params.UM1_I = struct( ...
    'scanner',           'Siemens MAGNETOM Allegra syngo MR A30 3T', ...
    'field_strength_T',  3.0, ...
    'TR',                1.500, ...   % confirmed: TR 1500 ms
    'TE_s',              0.025, ...   % confirmed: TE 25 ms
    'flip_angle_deg',    70, ...      % confirmed: 70 deg
    'nslices',           29, ...      % confirmed: Slices 29
    'slice_thickness_mm',4.0, ...     % confirmed: Slice thickness 4 mm
    'voxel_mm',          [3.1 3.1 4.0], ... % confirmed: 3.1×3.1×4.0 mm
    'slice_order',       'interleaved (Siemens)', ...  % confirmed: Multi-slice mode Interleaved; 29 slices odd-first
    'multiband_factor',  1, ...
    'nvolumes',          200, ...     % confirmed: Measurements 200
    'phase_enc_dir',     'A>>P', ...  % confirmed: Phase enc. dir. A >> P
    'SENSE_factor',      1 ...        % confirmed: PAT mode None
);

% ------------------------------------------------------------------
% UM2_I — University of Michigan, Ann Arbor MI, USA — Scanner 2 (ABIDE I)
% Scanner : GE 3T (reverse spiral sequence)
% Source  : func_UM_2.pdf, anat_UM_2.pdf
% GE reverse spiral (non-Cartesian EPI); slice timing correction is an
% approximation — 'ascending' used here as a conservative default.
% No SENSE/GRAPPA acceleration. 300 volumes × TR 2s = 10 min scan.
% Voxel 3.438×3.438×3.0 mm (FOV 220mm / base resolution 64 = 3.4375mm in-plane).
% ------------------------------------------------------------------
% site_params.UM2_I = struct( ...
%     'scanner',           'GE 3T (reverse spiral)', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.000, ...
%     'TE_s',              0.030, ...           % corrected from PDF (was 0.028)
%     'flip_angle_deg',    90, ...
%     'nslices',           40, ...              % corrected from PDF (was 38)
%     'slice_thickness_mm',3.0, ...             % corrected from PDF (was 3.5)
%     'voxel_mm',          [3.438 3.438 3.0], ...% corrected from PDF (was [3.0 3.0 3.5])
%     'slice_order',       'ascending', ...     % approx; GE reverse spiral non-standard
%     'multiband_factor',  1, ...
%     'nvolumes',          300, ...             % corrected from PDF (was 150)
%     'phase_enc_dir',     'A>>P', ...
%     'SENSE_factor',      1 ...                % GE scanner, no SENSE (was 2)
% );

% ------------------------------------------------------------------
% USM_I — University of Utah School of Medicine, Salt Lake City UT, USA (ABIDE I)
% Scanner : Siemens MAGNETOM TrioTim syngo MR B17 3T
% Source  : func_USM.pdf, anat_USM.pdf
% NOTE: Different from ABIDE II site_params.USM (which uses a placeholder).
%   ABIDE I USM: 40 slices, GRAPPA×2, TR=2s, TE=28ms — confirmed from PDF.
% 40 slices (even) → Siemens interleaved even-first.
% 240 volumes × TR 2s = 8 min scan.
% ------------------------------------------------------------------
% site_params.USM_I = struct( ...
%     'scanner',           'Siemens MAGNETOM TrioTim syngo MR B17 3T', ...
%     'field_strength_T',  3.0, ...
%     'TR',                2.000, ...
%     'TE_s',              0.028, ...
%     'flip_angle_deg',    90, ...
%     'nslices',           40, ...
%     'slice_thickness_mm',3.0, ...
%     'voxel_mm',          [3.4 3.4 3.0], ...
%     'slice_order',       'interleaved (Siemens)', ...   % 40 slices, even-first
%     'multiband_factor',  1, ...
%     'nvolumes',          240, ...
%     'phase_enc_dir',     'A>>P', ...
%     'SENSE_factor',      2 ...    % GRAPPA×2 (PAT mode GRAPPA, Accel. factor PE=2)
% );

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

% Enumerate all institute directories (ABIDE I and II)
inst_dirs = [dir(fullfile(imageFilesDir, 'ABIDEI-*')); ...
             dir(fullfile(imageFilesDir, 'ABIDEII-*'))];
inst_dirs = inst_dirs([inst_dirs.isdir]);


% to set starting subject according to previous stop point
start_from = 310;
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

            % Prefer rest_1; fall back to the first available run
            if ismember('rest_1', rest_names)
                selected_rest = 'rest_1';
            else
                selected_rest = rest_names{1};
                fprintf('    [WARN] rest_1 absent; using %s\n', selected_rest);
            end

            func_nii = find_nii(fullfile(ses_path, selected_rest), 'rest');
            if isempty(func_nii)
                warning('    No NIfTI in %s — skipping.', selected_rest);
                failed_subjects{end+1} = sprintf('%s_%s', subj_id, ses_name);
                continue;
            end
            FUNCTIONAL_FILES = {func_nii};
            nruns = 1;

            % ---------------------------------------------------------
            % CONN PROJECT FILE
            % ---------------------------------------------------------
            proj_name = sprintf('%s_%s.mat', subj_id, ses_name);
            proj_path = fullfile(siteOutDir, proj_name);

            if exist(proj_path, 'file')
                proj_dir = strrep(proj_path, '.mat', '');
                if exist(fullfile(proj_dir, 'results', 'firstlevel'), 'dir')
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
