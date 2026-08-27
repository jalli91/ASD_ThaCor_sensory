% ABIDE_motion_report.m
%
% Collects motion statistics for every preprocessed ABIDE subject and writes
% a summary Excel file.
%
% For each subject the script reads:
%   rp_*.txt                        — SPM realignment parameters (6 columns)
%   art_regression_outliers_*.mat   — ART scrubbing indicator matrix (R)
%
% Phenotypic CSV files (AdultFileID_ByGender/) are read to obtain the
% canonical SITE_ID string for each subject and to build one-hot institute
% columns in the output table.
%
% Computed metrics:
%   MeanFD          : mean framewise displacement (mm); rotations converted
%                     to mm using a 50-mm sphere radius, identical to the
%                     reference script ABIDE_motion_analx.m
%   ScrubFrames     : number of volumes flagged by ART
%   ScrubProportion : ScrubFrames / TotalFrames  (site-specific denominator
%                     taken from the row count of rp_*.txt)
%
% Folder naming convention (applied by rename_subjects.py / add_group.sh):
%   1st char  :  m = male (SEX 1)   |  f = female (SEX 2)
%   2nd char  :  a = ASD  (DX 1)    |  c = control (DX 2)
%   remaining :  numeric subject ID
%   Example   :  ma0050952 = male ASD subject 0050952
%
% Subjects with multiple sessions or rest runs are given one row per run;
% the SubjectID column always contains just the numeric subject ID.
%
% Output
%   /Volumes/Milk/CYL/ASD_ThaCor/ABIDE/motion_report.xlsx
%   Columns: Institute | SiteID | SubjectID | Sex | DX_GROUP |
%            Age | FIQ | Session | RestRun | MeanFD | ScrubFrames | ScrubProportion |
%            Site_<name1> | Site_<name2> | ...  (one-hot per institute)
%
%   Session  : numeric index extracted from the session_N folder name
%   RestRun  : numeric index extracted from the rest_N folder name

clear all;

%% ============================================================
%  CONFIGURATION
%% ============================================================
imageFilesDir = '/Volumes/Milk/CYL/ASD_ThaCor/ABIDE/ImageFiles';
phenoDir      = '/Volumes/Milk/CYL/ASD_ThaCor/ABIDE/AdultFileID_ByGender';
outputFile    = '/Volumes/Milk/CYL/ASD_ThaCor/ABIDE/profile_ABIDE.xlsx';

r = 50;   % sphere radius (mm) for converting rotations to mm displacement

%% ============================================================
%  LOAD PHENOTYPIC CSVs — build SubjectID → SITE_ID lookup
%% ============================================================
csv_files = {
    fullfile(phenoDir, 'ABIDEI_adult_male_FIQ70plus.csv');
    fullfile(phenoDir, 'ABIDEI_adult_female_FIQ70plus.csv');
    fullfile(phenoDir, 'ABIDEII_adult_male_FIQ70plus.csv');
    fullfile(phenoDir, 'ABIDEII_adult_female_FIQ70plus.csv')
};

% Maps: numeric SUB_ID (double) → phenotypic values
subj2site = containers.Map('KeyType', 'double', 'ValueType', 'char');
subj2age  = containers.Map('KeyType', 'double', 'ValueType', 'double');
subj2fiq  = containers.Map('KeyType', 'double', 'ValueType', 'double');

% SRS maps (raw scores; NaN when not available)
% ABIDE I  column names: SRS_RAW_TOTAL, SRS_AWARENESS, SRS_COGNITION,
%                        SRS_COMMUNICATION, SRS_MOTIVATION, SRS_MANNERISMS
% ABIDE II column names: SRS_TOTAL_RAW, SRS_AWARENESS_RAW, SRS_COGNITION_RAW,
%                        SRS_COMMUNICATION_RAW, SRS_MOTIVATION_RAW, SRS_MANNERISMS_RAW
SRS_FIELDS = {
    'SRS_Total',         {'SRS_RAW_TOTAL',    'SRS_TOTAL_RAW'};
    'SRS_Awareness',     {'SRS_AWARENESS',    'SRS_AWARENESS_RAW'};
    'SRS_Cognition',     {'SRS_COGNITION',    'SRS_COGNITION_RAW'};
    'SRS_Communication', {'SRS_COMMUNICATION','SRS_COMMUNICATION_RAW'};
    'SRS_Motivation',    {'SRS_MOTIVATION',   'SRS_MOTIVATION_RAW'};
    'SRS_Mannerisms',    {'SRS_MANNERISMS',   'SRS_MANNERISMS_RAW'};
};
n_srs = size(SRS_FIELDS, 1);
subj2srs = cell(n_srs, 1);
for k = 1:n_srs
    subj2srs{k} = containers.Map('KeyType', 'double', 'ValueType', 'double');
end

for f = 1:length(csv_files)
    if ~isfile(csv_files{f})
        fprintf('[WARN] Phenotypic CSV not found — skipping: %s\n', csv_files{f});
        continue;
    end

    % Read only the columns we need; tolerates extra/missing columns
    opts = detectImportOptions(csv_files{f}, 'Delimiter', ',');
    available = opts.VariableNames;
    available_upper = upper(strtrim(available));   % normalised for comparison

    % Required: SITE_ID, SUB_ID
    site_col_mask = strcmp(available_upper, 'SITE_ID');
    sub_col_mask  = strcmp(available_upper, 'SUB_ID');
    if ~any(site_col_mask) || ~any(sub_col_mask)
        fprintf('[WARN] SITE_ID or SUB_ID column not found in %s — skipping.\n', csv_files{f});
        continue;
    end

    % Optional but expected: AGE_AT_SCAN, FIQ
    age_col_mask = strcmp(available_upper, 'AGE_AT_SCAN');
    fiq_col_mask = strcmp(available_upper, 'FIQ');

    % Resolve SRS columns: pick the first candidate present for each subscale
    srs_col_found = cell(n_srs, 1);   % actual column name in this CSV (or '')
    for k = 1:n_srs
        candidates = SRS_FIELDS{k, 2};
        for c = 1:length(candidates)
            if any(strcmp(available_upper, candidates{c}))
                srs_col_found{k} = available{strcmp(available_upper, candidates{c})};
                break;
            end
        end
    end

    sel_cols = {available{site_col_mask}, available{sub_col_mask}};
    if any(age_col_mask), sel_cols{end+1} = available{age_col_mask}; end
    if any(fiq_col_mask), sel_cols{end+1} = available{fiq_col_mask}; end
    for k = 1:n_srs
        if ~isempty(srs_col_found{k}), sel_cols{end+1} = srs_col_found{k}; end
    end

    opts.SelectedVariableNames = sel_cols;
    T_pheno = readtable(csv_files{f}, opts);

    % Normalise column names to remove any trailing whitespace artefacts
    T_pheno.Properties.VariableNames = upper(strtrim(T_pheno.Properties.VariableNames));

    has_age = ismember('AGE_AT_SCAN', T_pheno.Properties.VariableNames);
    has_fiq = ismember('FIQ',         T_pheno.Properties.VariableNames);

    n_added = 0;
    for row = 1:height(T_pheno)
        raw_sub = T_pheno.SUB_ID(row);
        if iscell(raw_sub)
            sid = str2double(raw_sub{1});
        else
            sid = double(raw_sub);
        end

        if isnan(sid) || sid <= 0, continue; end

        raw_site = T_pheno.SITE_ID(row);
        if iscell(raw_site)
            site_str = strtrim(raw_site{1});
        else
            site_str = strtrim(char(raw_site));
        end

        if ~isempty(site_str)
            subj2site(sid) = site_str;
        end

        if has_age
            raw_age = T_pheno.AGE_AT_SCAN(row);
            if iscell(raw_age), raw_age = raw_age{1}; end
            age_val = double(raw_age);
            if ~isnan(age_val), subj2age(sid) = age_val; end
        end

        if has_fiq
            raw_fiq = T_pheno.FIQ(row);
            if iscell(raw_fiq), raw_fiq = raw_fiq{1}; end
            fiq_val = double(raw_fiq);
            if ~isnan(fiq_val) && fiq_val > 0, subj2fiq(sid) = fiq_val; end
        end

        for k = 1:n_srs
            if isempty(srs_col_found{k}), continue; end
            col_upper = upper(strtrim(srs_col_found{k}));
            if ~ismember(col_upper, T_pheno.Properties.VariableNames), continue; end
            raw_srs = T_pheno.(col_upper)(row);
            if iscell(raw_srs), raw_srs = raw_srs{1}; end
            if ischar(raw_srs) || isstring(raw_srs)
                srs_val = str2double(raw_srs);
            else
                srs_val = double(raw_srs);
            end
            if isscalar(srs_val) && ~isnan(srs_val) && srs_val >= 0
                subj2srs{k}(sid) = srs_val;
            end
        end

        n_added = n_added + 1;
    end
    fprintf('Loaded %4d subject entries (site/age/FIQ) from %s\n', n_added, csv_files{f});
end

srs_counts = cellfun(@(m) m.Count, subj2srs);
fprintf('Total subject→site mappings: %d  |  age: %d  |  FIQ: %d\n', ...
        subj2site.Count, subj2age.Count, subj2fiq.Count);
for k = 1:n_srs
    fprintf('  SRS %-20s: %d subjects\n', SRS_FIELDS{k,1}, srs_counts(k));
end
fprintf('\n');

%% ============================================================
%  ENUMERATE INSTITUTES
%% ============================================================
inst_dirs = [dir(fullfile(imageFilesDir, 'ABIDEI-*')); ...
             dir(fullfile(imageFilesDir, 'ABIDEII-*'))];
inst_dirs = inst_dirs([inst_dirs.isdir]);

if isempty(inst_dirs)
    error('No ABIDEI-* or ABIDEII-* directories found under:\n  %s', imageFilesDir);
end

% Pre-allocate result storage
result_institute    = {};
result_site_id      = {};   % canonical SITE_ID from phenotypic CSV
result_subjectID    = {};
result_sex          = [];
result_dx           = [];
result_age          = [];   % AGE_AT_SCAN from phenotypic CSV (NaN if missing)
result_fiq          = [];   % FIQ from phenotypic CSV (NaN if missing)
result_session      = [];   % numeric session number (from session_N folder name)
result_run          = [];   % numeric rest-run number (from rest_N folder name)
result_meanFD       = [];
result_scrubFrames  = [];
result_scrubProp    = [];
result_srs          = cell(n_srs, 1);   % one numeric array per SRS field
for k = 1:n_srs, result_srs{k} = []; end

total_processed = 0;
total_skipped   = 0;

%% ============================================================
%  MAIN LOOP
%% ============================================================
for inst_idx = 1:length(inst_dirs)

    inst_name = inst_dirs(inst_idx).name;   % e.g. 'ABIDEI-NYU'
    inst_path = fullfile(imageFilesDir, inst_name);

    fprintf('\n=== %s ===\n', inst_name);

    % Collect subject folders with recognised prefixes
    entries    = dir(inst_path);
    entries    = entries([entries.isdir]);
    all_names  = {entries.name};
    subj_mask  = (startsWith(all_names, 'ma') | startsWith(all_names, 'mc') | ...
                  startsWith(all_names, 'fa') | startsWith(all_names, 'fc')) & ...
                 ~contains(all_names, '_lowIQ');
    subj_names = all_names(subj_mask);

    if isempty(subj_names)
        fprintf('  [INFO] No renamed subject folders found (ma/mc/fa/fc) — skipping.\n');
        continue;
    end

    for subj_idx = 1:length(subj_names)

        subj_folder = subj_names{subj_idx};
        subj_path   = fullfile(inst_path, subj_folder);

        % ---- Decode folder name ----------------------------------------
        gender_char = subj_folder(1);   % 'm' | 'f'
        group_char  = subj_folder(2);   % 'a' (ASD) | 'c' (control)
        subj_id_num = subj_folder(3:end);   % e.g. '0050952'

        sex      = 1 + (gender_char == 'f');   % m→1, f→2
        dx_group = 1 + (group_char  == 'c');   % a→1, c→2

        % ---- Resolve phenotypic values from CSV ------------------------
        subj_num = str2double(subj_id_num);   % strips leading zeros

        if isKey(subj2site, subj_num)
            site_id = subj2site(subj_num);
        else
            site_id = inst_name;   % fall back to folder-level institute name
            fprintf('  [WARN] Subject %s not found in phenotypic CSVs — site set to folder name\n', ...
                    subj_id_num);
        end

        if isKey(subj2age, subj_num)
            subj_age = subj2age(subj_num);
        else
            subj_age = NaN;
        end

        if isKey(subj2fiq, subj_num)
            subj_fiq = subj2fiq(subj_num);
        else
            subj_fiq = NaN;
        end

        subj_srs = NaN(n_srs, 1);
        for k = 1:n_srs
            if isKey(subj2srs{k}, subj_num)
                subj_srs(k) = subj2srs{k}(subj_num);
            end
        end

        % ---- Discover sessions -----------------------------------------
        sess_dirs = dir(fullfile(subj_path, 'session_*'));
        sess_dirs = sess_dirs([sess_dirs.isdir]);

        if isempty(sess_dirs)
            fprintf('  [SKIP] %-20s no session_* folder\n', subj_folder);
            total_skipped = total_skipped + 1;
            continue;
        end

        for ses_idx = 1:length(sess_dirs)

            ses_name = sess_dirs(ses_idx).name;
            ses_path = fullfile(subj_path, ses_name);

            % ---- Discover rest runs ------------------------------------
            rest_dirs = dir(fullfile(ses_path, 'rest_*'));
            rest_dirs = rest_dirs([rest_dirs.isdir]);

            if isempty(rest_dirs)
                fprintf('  [SKIP] %-20s  %s: no rest_* folder\n', subj_folder, ses_name);
                total_skipped = total_skipped + 1;
                continue;
            end

            for run_idx = 1:length(rest_dirs)

                run_name = rest_dirs(run_idx).name;
                run_path = fullfile(ses_path, run_name);

                % ========================================================
                %  1. MOTION PARAMETERS  (rp_*.txt)
                % ========================================================
                rp_files = dir(fullfile(run_path, 'rp_*.txt'));

                if isempty(rp_files)
                    fprintf('  [SKIP] %-20s  %s/%s: rp_*.txt not found\n', ...
                            subj_folder, ses_name, run_name);
                    total_skipped = total_skipped + 1;
                    continue;
                end

                rp_path     = fullfile(rp_files(1).folder, rp_files(1).name);
                motion      = load(rp_path);     % N x 6  [x y z pitch roll yaw]
                total_frames = size(motion, 1);   % site-specific frame count

                % Framewise displacement (Power 2012 method)
                d_motion        = diff(motion, 1, 1);      % (N-1) x 6
                d_motion(:,4:6) = d_motion(:,4:6) * r;     % radians → mm
                FD              = sum(abs(d_motion), 2);    % (N-1) x 1
                meanFD          = mean(FD);

                % ========================================================
                %  2. SCRUBBING INDICATOR  (art_regression_outliers_au*.mat)
                % ========================================================
                % ART writes R as a T×K binary indicator matrix: one column
                % per flagged volume, with a single 1 at that timepoint.
                % The number of scrubbed frames is therefore size(R,2).
                % When no volumes are flagged ART may write an empty T×0
                % matrix, giving size(R,2)=0 (zero scrubbed frames).
                %
                % Search in run folder first, then session folder.
                art_files = dir(fullfile(run_path, 'art_regression_outliers_au*.mat'));
                if isempty(art_files)
                    art_files = dir(fullfile(ses_path, 'art_regression_outliers_au*.mat'));
                end

                if isempty(art_files)
                    fprintf('  [WARN] %-20s  %s/%s: art file not found — scrub=NaN\n', ...
                            subj_folder, ses_name, run_name);
                    scrub_frames = NaN;
                    scrub_prop   = NaN;
                else
                    art_path = fullfile(art_files(1).folder, art_files(1).name);
                    art_data = load(art_path);

                    % Field may be named 'R' (standard) or the first variable
                    if isfield(art_data, 'R')
                        R = art_data.R;
                    else
                        flds = fieldnames(art_data);
                        R    = art_data.(flds{1});
                    end

                    scrub_frames = size(R, 2);   % number of flagged-volume columns
                    scrub_prop   = scrub_frames / total_frames;
                end

                % ========================================================
                %  3. STORE RESULT
                % ========================================================
                % Extract numeric index from session_N / rest_N folder names
                ses_num_str = regexp(ses_name, '\d+', 'match', 'once');
                run_num_str = regexp(run_name, '\d+', 'match', 'once');
                ses_num = str2double(ses_num_str);   % NaN if no digits found
                run_num = str2double(run_num_str);

                result_institute{end+1,1}   = inst_name;
                result_site_id{end+1,1}     = site_id;
                result_subjectID{end+1,1}   = subj_id_num;
                result_sex(end+1,1)         = sex;
                result_dx(end+1,1)          = dx_group;
                result_age(end+1,1)         = subj_age;
                result_fiq(end+1,1)         = subj_fiq;
                result_session(end+1,1)     = ses_num;
                result_run(end+1,1)         = run_num;
                result_meanFD(end+1,1)      = meanFD;
                result_scrubFrames(end+1,1) = scrub_frames;
                result_scrubProp(end+1,1)   = scrub_prop;
                for k = 1:n_srs
                    result_srs{k}(end+1,1) = subj_srs(k);
                end

                fprintf('  %-22s  Site=%-18s  Sex=%d DX=%d  Age=%s  FIQ=%s  MeanFD=%.4f  Scrub=%s/%d  Prop=%s\n', ...
                        subj_id_num, site_id, sex, dx_group, ...
                        num_or_nan(subj_age), num_or_nan(subj_fiq), meanFD, ...
                        num_or_nan(scrub_frames), total_frames, ...
                        num_or_nan(scrub_prop));

                total_processed = total_processed + 1;

            end   % run loop
        end   % session loop
    end   % subject loop
end   % institute loop

%% ============================================================
%  WRITE EXCEL REPORT
%% ============================================================
fprintf('\n--- Summary: %d rows processed, %d subjects/runs skipped ---\n', ...
        total_processed, total_skipped);

if total_processed == 0
    fprintf('[WARNING] No preprocessed data found. Excel file not written.\n');
    fprintf('Make sure preprocessing has been run and subject folders are renamed\n');
    fprintf('(ma/mc/fa/fc prefix) before running this script.\n');
    return;
end

%% ============================================================
%  BUILD DUMMY-CODED INSTITUTE VARIABLES  (reference = ABIDEII-OILH_2)
%% ============================================================
BASELINE_SITE = 'ABIDEII-OILH_2';   % reference category — all dummies = 0

% Collect all unique SITE_IDs observed in this run (alphabetical order)
unique_sites = unique(result_site_id);
unique_sites = sort(unique_sites);
n_sites = length(unique_sites);
n_rows  = length(result_site_id);

% Verify the baseline exists in the data
baseline_idx = find(strcmp(unique_sites, BASELINE_SITE));
if isempty(baseline_idx)
    fprintf('[WARN] Baseline site "%s" not found in data — falling back to full dummy set.\n', ...
            BASELINE_SITE);
    dummy_sites = unique_sites;          % include all (no true baseline)
else
    fprintf('\nBaseline (reference) site: %s  (%d subjects/runs)\n', ...
            BASELINE_SITE, sum(strcmp(result_site_id, BASELINE_SITE)));
    % Exclude the baseline; remaining sites each get one dummy column
    dummy_sites = unique_sites([1:baseline_idx-1, baseline_idx+1:end]);
end

n_dummy = length(dummy_sites);
n_sites_display = n_sites;   % keep for summary print

fprintf('\n--- Institute dummy coding: %d sites total, %d dummy columns (baseline excluded) ---\n', ...
        n_sites_display, n_dummy);
for s = 1:n_sites
    n_in_site = sum(strcmp(result_site_id, unique_sites{s}));
    is_base   = strcmp(unique_sites{s}, BASELINE_SITE);
    fprintf('  %-25s  %d rows%s\n', unique_sites{s}, n_in_site, ...
            repmat('  [BASELINE]', 1, is_base));
end

% Build binary matrix  (n_rows × n_dummy)  — baseline rows are all-zero
onehot_matrix = zeros(n_rows, n_dummy);
for s = 1:n_dummy
    onehot_matrix(:, s) = strcmp(result_site_id, dummy_sites{s});
end

% Derive valid MATLAB table variable names:  'Site_CALTECH', 'Site_ABIDEII_BNI_1', …
% makeValidName can map two distinct SITE_IDs to the same string (e.g. 'UM-1' and
% 'UM_1' both become 'UM_1').  makeUniqueStrings appends _2, _3, … to resolve any
% collisions so every dummy column has a distinct name in the output table.
site_varnames = cellfun(@(s) ['Site_' matlab.lang.makeValidName(s)], ...
                        dummy_sites, 'UniformOutput', false);
site_varnames = matlab.lang.makeUniqueStrings(site_varnames);

%% ============================================================
%  ASSEMBLE TABLE
%% ============================================================
T = table(result_site_id, result_subjectID, ...
          result_sex, result_dx, result_age, result_fiq, ...
          result_session, result_run, ...
          result_meanFD, result_scrubFrames, result_scrubProp, ...
          'VariableNames', {'Site', 'ID', 'Sex', 'GROUP', ...
                            'Age', 'FSIQ', 'Session', 'RestRun', ...
                            'MeanFD', 'ScrubFrames', 'ScrubProportion'});

% Append SRS columns (Total first, then subscores in canonical order)
for k = 1:n_srs
    T.(SRS_FIELDS{k,1}) = result_srs{k};
end

% Flag subjects with ScrubProportion >= 20% (placed right after ScrubProportion)
SCRUB_THRESHOLD = 0.20;
high_motion_flag = repmat({''}, height(T), 1);
high_motion_flag(~isnan(T.ScrubProportion) & T.ScrubProportion >= SCRUB_THRESHOLD) = {'x'};
T.HighMotion = high_motion_flag;

% Flag subjects whose MeanFD exceeds mean + 2*SD across all subjects
valid_fd      = T.MeanFD(~isnan(T.MeanFD));
fd_mean       = mean(valid_fd);
fd_std        = std(valid_fd);
fd_threshold  = fd_mean + 2 * fd_std;

high_fd_flag  = repmat({''}, height(T), 1);
high_fd_flag(~isnan(T.MeanFD) & T.MeanFD > fd_threshold) = {'x'};
T.HighMeanFD  = high_fd_flag;

fprintf('MeanFD stats (all subjects): mean=%.4f  SD=%.4f  threshold(mean+2SD)=%.4f\n', ...
        fd_mean, fd_std, fd_threshold);
fprintf('Subjects exceeding MeanFD threshold: %d / %d\n', ...
        sum(strcmp(T.HighMeanFD, 'x')), height(T));

% Append dummy columns (baseline site has no column — all-zero by convention)
for s = 1:n_dummy
    T.(site_varnames{s}) = onehot_matrix(:, s);
end

% Sort by institute, subject ID, session, then run for readability
T = sortrows(T, {'Site', 'ID', 'Session', 'RestRun'});

writetable(T, outputFile, 'Sheet', 'MotionReport', 'WriteMode', 'overwritesheet');
fprintf('\nReport written to: %s\n', outputFile);
srs_col_list = strjoin(cellfun(@(f) f, SRS_FIELDS(:,1), 'UniformOutput', false), ' | ');
fprintf('Fixed columns : Site | ID | Sex (1=M,2=F) | GROUP (1=ASD,2=TDC) | Age | FIQ | Session | RestRun | MeanFD | ScrubFrames | ScrubProportion | HighMotion | HighMeanFD | %s\n', srs_col_list);
fprintf('One-hot cols  : %s\n', strjoin(site_varnames, ' | '));

%% ============================================================
%  HIGHLIGHT HIGH-MOTION ROWS (ScrubProportion >= 20%)
%% ============================================================
flag_mask  = strcmp(T.HighMotion, 'x');
flag_rows  = find(flag_mask);   % 1-based row indices in the sorted table
excel_rows = flag_rows + 1;     % +1 to skip the Excel header row

fprintf('\nSubjects with ScrubProportion >= 20%%: %d\n', numel(excel_rows));

if ~isempty(excel_rows)
    % Build a Python list literal: [2, 5, 10, ...]
    row_list_str = ['[' strjoin(arrayfun(@num2str, excel_rows, 'UniformOutput', false), ', ') ']'];

    python_code = sprintf([ ...
        'import openpyxl\n' ...
        'from openpyxl.styles import PatternFill\n' ...
        'wb = openpyxl.load_workbook(r"%s")\n' ...
        'ws = wb["MotionReport"]\n' ...
        'fill = PatternFill(start_color="FFFF00", end_color="FFFF00", fill_type="solid")\n' ...
        'for r in %s:\n' ...
        '    for c in range(1, ws.max_column + 1):\n' ...
        '        ws.cell(row=r, column=c).fill = fill\n' ...
        'wb.save(r"%s")\n' ...
    ], outputFile, row_list_str, outputFile);

    try
        pyrun(python_code);
        fprintf('Yellow row highlighting applied to %d subjects (ScrubProportion >= 20%%).\n', numel(excel_rows));
    catch ME
        fprintf('[WARN] Could not apply Excel row highlighting: %s\n', ME.message);
        fprintf('       Ensure Python is configured in MATLAB and openpyxl is installed.\n');
        fprintf('       To install: pip install openpyxl\n');
    end
end

%% ============================================================
%  HELPER
%% ============================================================
function s = num_or_nan(v)
% Return a display string — 'NaN' if v is NaN, otherwise the number.
    if isnan(v)
        s = 'NaN';
    else
        s = num2str(v);
    end
end
