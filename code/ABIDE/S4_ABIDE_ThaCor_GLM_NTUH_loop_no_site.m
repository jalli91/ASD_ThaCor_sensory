%  Loop over all profile_ABIDE+NTUH_matched*.xlsx files and run the full
%  thalamo-cortical GLM pipeline for each matched sample.
%
%  The CONN batch_output map is built ONCE before the loop (filesystem scan
%  independent of which profile is selected).  Per-profile work starts at
%  FC loading and ends at Excel export.
clear all
dataDir     = '/Volumes/Milk/CYL/ASD_ThaCor/ABIDE/ImageFiles';
profile_dir = '/Volumes/Milk/CYL/ASD_ThaCor/ABIDE/Processed_matched_files';
output_dir  = '/Volumes/Milk/CYL/ASD_ThaCor/ABIDE/GLM_results';

%% Find all matched profile Excel files
profile_hits = dir(fullfile(profile_dir, 'profile_ABIDE+fullNTUH_matched*.xlsx'));
if isempty(profile_hits)
    error('No profile_ABIDE+NTUH_matched*.xlsx found in %s', profile_dir);
end
fprintf('Found %d profile file(s) to process.\n\n', numel(profile_hits));

%% Build subject -> CONN results map ONCE (covers all available subjects)
% Subject folders are named [gender(1)][group(1)][numericID] (e.g. 'ma0050952')
% under each ABIDE*/ site directory; CONN outputs live in <site>/batch_output/.
inst_dirs = [dir(fullfile(dataDir, 'ABIDEI-*')); ...
             dir(fullfile(dataDir, 'ABIDEII-*')); ...
             dir(fullfile(dataDir, 'NTU*'))];
inst_dirs = inst_dirs([inst_dirs.isdir]);

batch_dirs = {};
for k = 1:length(inst_dirs)
    bd = fullfile(dataDir, inst_dirs(k).name, 'batch_output');
    if exist(bd, 'dir'), batch_dirs{end+1} = bd; end
end
ntuh_batch_dir = fullfile(dataDir, 'NTUH', 'batch_output');
if exist(ntuh_batch_dir, 'dir') && ~any(strcmp(batch_dirs, ntuh_batch_dir))
    batch_dirs{end+1} = ntuh_batch_dir;
end

subj_results_map = containers.Map('KeyType', 'double', 'ValueType', 'char');
for k = 1:length(batch_dirs)
    batch_dir = batch_dirs{k};
    entries   = dir(batch_dir);
    entries   = entries([entries.isdir]);
    for e = 1:length(entries)
        sname = entries(e).name;
        % Match folders named {prefix}{id}_session_{n}, e.g. ma2002_session_1
        tok = regexp(sname, '^(ma|fa|mc|fc)(\d+)_session_(\d+)$', 'tokens', 'once');
        if isempty(tok), continue; end
        numeric_id = str2double(tok{2});
        if isnan(numeric_id), continue; end
        sess_num = str2double(tok{3});
        if isKey(subj_results_map, numeric_id) && sess_num > 1, continue; end
        candidate = fullfile(batch_dir, sname, ...
            'results', 'firstlevel', 'SBC_01', ...
            'resultsROI_Subject001_Condition001.mat');
        if exist(candidate, 'file')
            subj_results_map(numeric_id) = candidate;
        end
    end
end
fprintf('CONN results map built: %d subjects indexed.\n\n', subj_results_map.Count);

%% ========== Main loop over profile files ==========
for pf = 1:numel(profile_hits)

fprintf('========== Profile %d/%d: %s ==========\n', pf, numel(profile_hits), profile_hits(pf).name);

profile_file = fullfile(profile_dir, profile_hits(pf).name);
[~, fname, ~] = fileparts(profile_hits(pf).name);          % e.g. profile_ABIDE+NTUH_matched_AM_3_0.2_y
match_suffix  = fname(length('use_profile_ABIDE+fullNTUH_matched')+1:end);  % e.g. _AM_3_0.2_y

%% Read matched subject profile from Excel
% Sheet 'matched_all': Sex = Male 1 / Female 0; GROUP = ASD 1 / Control 0
Profile_All = readtable(profile_file, 'Sheet', 'matched_all');
SubNum = height(Profile_All);
fprintf('  Subjects in profile: %d\n', SubNum);

%% Build covariate matrix: [age, meanFD, sex, group, sex*group, site dummies (all sites)]
% Columns 6-16: Site_ABIDEII_IP_1, Site_ABIDEII_IU_1, Site_ABIDEII_NYU_1, Site_ABIDEII_OILH_2,
%               Site_ABIDEII_USM_1, Site_CALTECH, Site_CMU, Site_MAX_MUN, Site_NYU, Site_OLIN, Site_UM_1
Profile_All_ext = zeros(SubNum, 19);
for i = 1:SubNum
    Profile_All_ext(i,  1) = Profile_All.Age(i);
    Profile_All_ext(i,  2) = Profile_All.MeanFD(i);
    Profile_All_ext(i,  3) = Profile_All.Sex(i);        % Male=1, Female=0
    Profile_All_ext(i,  4) = Profile_All.GROUP(i);      % ASD=1, Control=0
    Profile_All_ext(i,  5) = Profile_All_ext(i, 4) * Profile_All_ext(i, 3);  % sex*group
end

X = [ones(SubNum,1), Profile_All_ext];  % 加上 intercept
fprintf('Matrix size: %d x %d\n', size(X,1), size(X,2));
fprintf('Rank: %d\n', rank(X));

ntuh = all(Profile_All_ext(:,6:16) == 0, 2);
sex = Profile_All_ext(:,3);
grp = Profile_All_ext(:,4);
fprintf('NTUH Male-ASD: %d\n',   sum(ntuh & sex==1 & grp==1));
fprintf('NTUH Male-Ctrl: %d\n',  sum(ntuh & sex==1 & grp==0));
fprintf('NTUH Female-ASD: %d\n', sum(ntuh & sex==0 & grp==1));
fprintf('NTUH Female-Ctrl: %d\n',sum(ntuh & sex==0 & grp==0));

%% Load FC data for each matched subject (ordered by profile)
profile_ids = cellfun(@str2double, Profile_All.ID);
missing_ids = [];
FC_raw      = {};
names       = {};   % ROI names — captured from first successfully loaded .mat

for i = 1:SubNum
    numeric_id = profile_ids(i);
    if ~isKey(subj_results_map, numeric_id)
        missing_ids(end+1) = numeric_id; %#ok<AGROW>
        FC_raw{end+1} = [];   % placeholder to keep indexing aligned with valid_mask
        continue;
    end
    ToOpen = subj_results_map(numeric_id);
    disp(ToOpen);
    loaded_vars = load(ToOpen);   % loads Z (and names, etc.) without polluting workspace
    FC_raw{end+1} = loaded_vars.Z;
    if isempty(names) && isfield(loaded_vars, 'names')
        names = loaded_vars.names;
    end
end

if ~isempty(missing_ids)
    fprintf('\n  [WARNING] No CONN result for %d subject(s): ', numel(missing_ids));
    fprintf('%d ', missing_ids);
    fprintf('\n');
end

% Determine max grid size; pad all matrices to uniform size with NaN
loaded_only = FC_raw(~cellfun(@isempty, FC_raw));
max_dim1 = max(cellfun(@(x) size(x,1), loaded_only));
max_dim2 = max(cellfun(@(x) size(x,2), loaded_only));
if max_dim1 ~= min(cellfun(@(x) size(x,1), loaded_only)) || ...
   max_dim2 ~= min(cellfun(@(x) size(x,2), loaded_only))
    fprintf('  [WARNING] FC matrix size mismatch — padding smaller grids with NaN.\n');
end
for i = 1:length(FC_raw)
    if isempty(FC_raw{i}), continue; end
    [d1, d2] = size(FC_raw{i});
    if d1 < max_dim1 || d2 < max_dim2
        Z_pad = NaN(max_dim1, max_dim2);
        Z_pad(1:d1, 1:d2) = FC_raw{i};
        FC_raw{i} = Z_pad;
    end
end

valid_raw = FC_raw(~cellfun(@isempty, FC_raw));
FC_All    = cat(3, valid_raw{:});

%% GLM: FC ~ age + meanFD + sex + group + sex*group + site_dummies
valid_mask = true(SubNum, 1);
for mi = 1:length(missing_ids)
    valid_mask(profile_ids == missing_ids(mi)) = false;
end
valid_idx    = find(valid_mask);
SubNum_valid = length(valid_idx);
[dim1, dim2] = deal(max_dim1, max_dim2);

FC_test          = {};
FC_test_sex      = {};
FC_test_interact = {};   % SexxGroup

for i = 1:dim1
    for j = 1:dim2
        if i ~= j
            FC_in = reshape(FC_All(i,j,:), [SubNum_valid, 1]);
            par_in = fitlm(Profile_All_ext(valid_idx,:), FC_in);
            interval = coefCI(par_in, 0.05);
            interval_value = par_in.Coefficients.Estimate(:) - interval(:,1);
            %[beta interval p-value] for group
            FC_test{i,j}          = [par_in.Coefficients.Estimate(5) interval_value(5) par_in.Coefficients.pValue(5)];
            %[beta interval p-value] for sex
            FC_test_sex{i,j}      = [par_in.Coefficients.Estimate(4) interval_value(4) par_in.Coefficients.pValue(4)];
            %[beta interval p-value] for sex*group interaction
            FC_test_interact{i,j} = [par_in.Coefficients.Estimate(6) interval_value(6) par_in.Coefficients.pValue(6)];
        else
            FC_test{i,j}          = nan;
            FC_test_sex{i,j}      = nan;
            FC_test_interact{i,j} = nan;
        end
    end
end

%% Within-sex subgroup GLM: FC ~ age + meanFD + group + site_dummies
% Split valid subjects into male (Sex=1) and female (Sex=0)
sex_valid       = Profile_All_ext(valid_idx, 3);
male_in_valid   = find(sex_valid == 1);
female_in_valid = find(sex_valid == 0);

% Reduced covariate matrix: drop sex (col 3) and sex*group (col 5)
% Columns: [age(1), meanFD(2), group(3), site_dummies(4-13)]
% fitlm coefficients: intercept(1), age(2), meanFD(3), group(4), site_dummies(5-14)
Profile_nosex      = [Profile_All_ext(valid_idx, 1:2), ...
                      Profile_All_ext(valid_idx, 4)];
Profile_male_ext   = Profile_nosex(male_in_valid,   :);
Profile_female_ext = Profile_nosex(female_in_valid, :);

FC_male   = FC_All(:, :, male_in_valid);
FC_female = FC_All(:, :, female_in_valid);

SubNum_male   = length(male_in_valid);
SubNum_female = length(female_in_valid);
fprintf('  Within-sex subgroups — Male: %d, Female: %d\n', SubNum_male, SubNum_female);

FC_test_male   = {};
FC_test_female = {};
for i = 1:dim1
    for j = 1:dim2
        if i ~= j
            % Male subgroup
            FC_in_m = reshape(FC_male(i,j,:),   [SubNum_male,   1]);
            par_m   = fitlm(Profile_male_ext,   FC_in_m);
            int_m   = coefCI(par_m, 0.05);
            iv_m    = par_m.Coefficients.Estimate(:) - int_m(:,1);
            % group is coefficient 4 in the reduced model
            FC_test_male{i,j} = [par_m.Coefficients.Estimate(4) iv_m(4) par_m.Coefficients.pValue(4)];

            % Female subgroup
            FC_in_f = reshape(FC_female(i,j,:), [SubNum_female, 1]);
            par_f   = fitlm(Profile_female_ext, FC_in_f);
            int_f   = coefCI(par_f, 0.05);
            iv_f    = par_f.Coefficients.Estimate(:) - int_f(:,1);
            FC_test_female{i,j} = [par_f.Coefficients.Estimate(4) iv_f(4) par_f.Coefficients.pValue(4)];
        else
            FC_test_male{i,j}   = nan;
            FC_test_female{i,j} = nan;
        end
    end
end

%% ROI index / name extraction (using names from loaded .mat)
% Somatosensory
ROI_soma_index = [];
ROI_soma_name  = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'PostCG')))        || ...
       ~isempty(cell2mat(strfind(names(i), 'Parietal Oper'))) || ...
       ~isempty(cell2mat(strfind(names(i), 'SPL')))           || ...
       ~isempty(cell2mat(strfind(names(i), 'aSMG')))          || ...
       ~isempty(cell2mat(strfind(names(i), 'pSMG')))
        ROI_soma_index = [ROI_soma_index, i];
        soma_name = char(names(i));
        ROI_soma_name = [ROI_soma_name, string(soma_name(15:end))];
    end
end

% Thalamus
ROI_tha_index = [];
ROI_tha_name  = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'Thalamus')))
        ROI_tha_index = [ROI_tha_index, i];
        thalamus_name = char(names(i));
        ROI_tha_name  = [ROI_tha_name, string(thalamus_name(15:end))];
    end
end
ROI_tha_name  = ["Thalamus R" "Thalamus L"];
ROI_soma_name = ["PostCG R" "PostCG L" "SPL R" "SPL L" "aSMG R" "aSMG L" "pSMG R" "pSMG L" "PO R" "PO L"];

% Auditory
ROI_au_index = [];
ROI_au_name  = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'HG')))   || ...
       ~isempty(cell2mat(strfind(names(i), 'pSTG'))) || ...
       ~isempty(cell2mat(strfind(names(i), 'aSTG')))
        Au_name     = char(names(i));
        Au_name_st  = cell2mat(strfind(names(i), '.'));
        Au_name_end = cell2mat(strfind(names(i), '('));
        Au_name_r   = cell2mat(strfind(names(i), 'r'));
        Au_name_l   = cell2mat(strfind(names(i), 'l'));
        if length(Au_name_st) == 1
            if Au_name_r(1) < Au_name_end
                ROI_au_index = [ROI_au_index, i];
                ROI_au_name  = [ROI_au_name, string([Au_name(Au_name_st+1:Au_name_r(1)-1) 'R'])];
            elseif Au_name_l(1) < Au_name_end
                ROI_au_index = [ROI_au_index, i];
                ROI_au_name  = [ROI_au_name, string([Au_name(Au_name_st+1:Au_name_l(3)-1) 'L'])];
            else
                ROI_au_index = [ROI_au_index, i];
                ROI_au_name  = [ROI_au_name, string(Au_name(Au_name_st+1:Au_name_end-1))];
            end
        end
    end
end

ROI_au_name = ["aSTG R", "aSTG L", "pSTG R", "pSTG L", "HG R", "HG L"];

% Visual
ROI_vi_index = [];
ROI_vi_name  = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'Visual')))
        Vi_name_st  = cell2mat(strfind(names(i), '.'));
        Vi_name_end = cell2mat(strfind(names(i), '('));
        if length(Vi_name_st) == 2
            Vi_name = char(names(i));
            ROI_vi_index = [ROI_vi_index, i];
            if length(Vi_name_end) == 2
                ROI_vi_name = [ROI_vi_name, string(Vi_name(Vi_name_st(1)+1:Vi_name_end(2)-1))];
            else
                ROI_vi_name = [ROI_vi_name, string(Vi_name(Vi_name_st(1)+1:Vi_name_end-1))];
            end
        end
    end
end

% Insular
ROI_ins_index = [];
ROI_ins_name  = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'Insular')))
        Ins_name_st  = cell2mat(strfind(names(i), '.'));
        Ins_name_end = cell2mat(strfind(names(i), '('));
        Ins_name = char(names(i));
        ROI_ins_index = [ROI_ins_index, i];
        ROI_ins_name  = [ROI_ins_name, string(Ins_name(Ins_name_st+1:Ins_name_end-1))];
    end
end
ROI_ins_name = ["IC R" "IC L"];

%% sub_extract for all contrasts
Test_ThaSoma          = sub_extract('ThaSoma group',        ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test);
Test_sex_ThaSoma      = sub_extract('ThaSoma sex',          ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_sex);
Test_interact_ThaSoma = sub_extract('ThaSoma group*sex',    ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_interact);
Test_male_ThaSoma     = sub_extract('ThaSoma male group',   ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_male);
Test_female_ThaSoma   = sub_extract('ThaSoma female group', ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_female);

Test_ThaAu            = sub_extract('ThaAu group',          ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test);
Test_sex_ThaAu        = sub_extract('ThaAu sex',            ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_sex);
Test_interact_ThaAu   = sub_extract('ThaAu group*sex',      ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_interact);
Test_male_ThaAu       = sub_extract('ThaAu male group',     ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_male);
Test_female_ThaAu     = sub_extract('ThaAu female group',   ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_female);

Test_ThaVi            = sub_extract('ThaVi group',          ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test);
Test_sex_ThaVi        = sub_extract('ThaVi sex',            ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_sex);
Test_interact_ThaVi   = sub_extract('ThaVi group*sex',      ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_interact);
Test_male_ThaVi       = sub_extract('ThaVi male group',     ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_male);
Test_female_ThaVi     = sub_extract('ThaVi female group',   ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_female);

Test_ThaIns           = sub_extract('ThaIns group',         ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test);
Test_sex_ThaIns       = sub_extract('ThaIns sex',           ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_sex);
Test_interact_ThaIns  = sub_extract('ThaIns group*sex',     ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_interact);
Test_male_ThaIns      = sub_extract('ThaIns male group',    ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_male);
Test_female_ThaIns    = sub_extract('ThaIns female group',  ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_female);

%% Export all Test_ results to Excel
fprintf('\n  Exporting results to Excel...\n');

outputFile = fullfile(output_dir, sprintf('ABIDE+fullNTUH_ThaCor_nosite%s.xlsx', match_suffix));
if exist(outputFile, 'file'), delete(outputFile); end

% Sheet definitions: each row is {sheet_name, {label, table; ...}}
sheet_defs = { ...
    'All', { ...
        'Test_ThaSoma',          Test_ThaSoma; ...
        'Test_sex_ThaSoma',      Test_sex_ThaSoma; ...
        'Test_interact_ThaSoma', Test_interact_ThaSoma; ...
        'Test_ThaAu',            Test_ThaAu; ...
        'Test_sex_ThaAu',        Test_sex_ThaAu; ...
        'Test_interact_ThaAu',   Test_interact_ThaAu; ...
        'Test_ThaVi',            Test_ThaVi; ...
        'Test_sex_ThaVi',        Test_sex_ThaVi; ...
        'Test_interact_ThaVi',   Test_interact_ThaVi; ...
        'Test_ThaIns',            Test_ThaIns; ...
        'Test_sex_ThaIns',       Test_sex_ThaIns; ...
        'Test_interact_ThaIns',  Test_interact_ThaIns ...
    }; ...
    'Female', { ...
        'Test_female_ThaSoma',   Test_female_ThaSoma; ...
        'Test_female_ThaAu',     Test_female_ThaAu; ...
        'Test_female_ThaVi',     Test_female_ThaVi; ...
        'Test_female_ThaIns',    Test_female_ThaIns ...
    }; ...
    'Male', { ...
        'Test_male_ThaSoma',     Test_male_ThaSoma; ...
        'Test_male_ThaAu',       Test_male_ThaAu; ...
        'Test_male_ThaVi',       Test_male_ThaVi; ...
        'Test_male_ThaIns',      Test_male_ThaIns ...
    } ...
};

for s = 1:size(sheet_defs, 1)
    sheet_name = sheet_defs{s, 1};
    tables     = sheet_defs{s, 2};
    combined   = {};
    cur_row    = 1;

    for t = 1:size(tables, 1)
        tbl_label = tables{t, 1};
        tbl       = tables{t, 2};
        [nr, nc]  = size(tbl);

        % Label row (separator between sub-tables)
        combined{cur_row, 1} = tbl_label;
        cur_row = cur_row + 1;

        % Write table contents
        for r = 1:nr
            for c = 1:nc
                v = tbl{r, c};
                if isnumeric(v) && isscalar(v) && ~isnan(v)
                    combined{cur_row + r - 1, c} = v;
                elseif isstring(v) || ischar(v)
                    combined{cur_row + r - 1, c} = char(v);
                else
                    combined{cur_row + r - 1, c} = '';
                end
            end
        end

        cur_row = cur_row + nr + 1;  % blank separator row between sub-tables
    end

    writecell(combined, outputFile, 'Sheet', sheet_name);
    fprintf('    Written sheet: %s\n', sheet_name);
end

fprintf('  Excel export complete: %s\n\n', outputFile);

end  % end profile loop

fprintf('All %d profile(s) processed.\n', numel(profile_hits));
