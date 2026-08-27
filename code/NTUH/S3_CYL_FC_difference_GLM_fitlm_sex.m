%  grouping FC by ROI*ROI into 3D array
clear all
dataDir = ['/Volumes/Milk/CYL/batch_output_0.01_0.08/'];

% Define the file path of gender according indexing file
filePath = '/Volumes/Milk/CYL/Profile_sheet/gender_index.xlsx';
% Read the Excel file into a table
gender_table = readtable(filePath);
gender_table = table2cell(gender_table);

%% grouping the result folder in the directory
file_list = dir(dataDir);
files_NAC = {};
files_ASD = {};
for i = 1:length(file_list)
    file_name = file_list(i).name;
    if startsWith(file_name, 'conn_1') & ~endsWith(file_name, '.mat')
        % NAC group file index start with 1
        files_NAC{end+1} = file_name;
    elseif startsWith(file_name, 'conn_2') & ~endsWith(file_name, '.mat')
        % ASD group file index start with 2
        files_ASD{end+1} = file_name;
    end
end
SubNum = length(files_NAC)+length(files_ASD);

% grouping file names by gender
files_NAC_m = {}; % collecting male NAC file
files_NAC_f = {}; % collecting female NAC file
for i = 1:length(files_NAC)
    for j = 1:length(gender_table)
        if sum(gender_table{j, 1}(2:5) == files_NAC{i}(6:9)) == 4
            if gender_table{j, 2} == 1 % M1F0
                files_NAC_m{end+1} = files_NAC{i};
            else
                files_NAC_f{end+1} = files_NAC{i};
            end
        end
    end
end    

files_ASD_m = {}; % collecting male ASD file
files_ASD_f = {}; % collecting female ASD file
for i = 1:length(files_ASD)
    for j = 1:length(gender_table)
        if sum(gender_table{j, 1}(2:5) == files_ASD{i}(6:9)) == 4
            if gender_table{j, 2} == 1 % M1F0
                files_ASD_m{end+1} = files_ASD{i};
            else
                files_ASD_f{end+1} = files_ASD{i};
            end
        end
    end
end 
SubNum = length(files_NAC)+length(files_ASD);
SubNum_m = length(files_NAC_m)+length(files_ASD_m);
SubNum_f = length(files_NAC_f)+length(files_ASD_f);

%% read Profile of age, meanFD and gender for the covariates matrix
Profile_All = readcell('/Volumes/Milk/CYL/Profile_sheet/All_profile.xlsx');
% extract age and meanFD only
% group by gender
n_m = 1; % male counter
n_f = 1; % female counter
for i = 1:SubNum
    if  Profile_All{i+1, 4} == 1
        % age
        Profile_All_ext_m(n_m, 1) = str2double(Profile_All{i+1, 2}(1:5));
        % meanFD
        Profile_All_ext_m(n_m, 2) = Profile_All{i+1, 3};
        % Group
        Profile_All_ext_m(n_m, 3) = Profile_All{i+1, 5};
        n_m = n_m + 1;
    elseif Profile_All{i+1, 4} == 0
        % age
        Profile_All_ext_f(n_f, 1) = str2double(Profile_All{i+1, 2}(1:5));
        % meanFD
        Profile_All_ext_f(n_f, 2) = Profile_All{i+1, 3};
        % Group
        Profile_All_ext_f(n_f, 3) = Profile_All{i+1, 5};
        n_f = n_f + 1;
    end
end

%% assign Z matrix from output into 3D array by group
% male
FC_NAC_m = [];
% for i = 1:length(files_NAC_m)
for i = 1:length(files_NAC_m)
    data_name = string(files_NAC_m(i));
    disp(data_name);
    ToOpen = append(dataDir, data_name, '/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat');
    load(ToOpen);
    FC_NAC_m = cat(3, FC_NAC_m, Z);
end      

FC_ASD_m = [];
% for i = 1:length(files_ASD_m)
for i = 1:length(files_ASD_m)
    data_name = string(files_ASD_m(i));
    disp(data_name);
    ToOpen = append(dataDir, data_name, '/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat');
    load(ToOpen);
    FC_ASD_m = cat(3, FC_ASD_m, Z);
end 
FC_All_m = cat(3, FC_NAC_m, FC_ASD_m);

% female
FC_NAC_f = [];
% for i = 1:length(files_NAC_f)
for i = 1:length(files_NAC_f)
    data_name = string(files_NAC_f(i));
    disp(data_name);
    ToOpen = append(dataDir, data_name, '/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat');
    load(ToOpen);
    FC_NAC_f = cat(3, FC_NAC_f, Z);
end  

FC_ASD_f = [];
% for i = 1:length(files_ASD_f)
for i = 1:length(files_ASD_f)
    data_name = string(files_ASD_f(i));
    disp(data_name);
    ToOpen = append(dataDir, data_name, '/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat');
    load(ToOpen);
    FC_ASD_f = cat(3, FC_ASD_f, Z);
end
FC_All_f = cat(3, FC_NAC_f, FC_ASD_f);

%% GLM with age, meanFD, gender and Group to test beta 
% male
FC_test_m = {};
[dim1, dim2] = size(Z);
for i = 1:dim1
    for j = 1:dim2
        if i ~= j
            FC_in = reshape(FC_All_m(i,j, :), [SubNum_m, 1]);
            par_in = fitlm(Profile_All_ext_m, FC_in);
            interval = coefCI(par_in, 0.05);
            interval_value = par_in.Coefficients.Estimate(:) - interval(:,1);
            %[beta p-value] for group
            FC_test_m{i, j} = [par_in.Coefficients.Estimate(4) interval_value(4) par_in.Coefficients.pValue(4)]; 
        else
            FC_test_m{i, j} = nan;
        end
    end    
end 

% female
FC_test_f = {};
[dim1, dim2] = size(Z);
for i = 1:dim1
    for j = 1:dim2
        if i ~= j
            FC_in = reshape(FC_All_f(i,j, :), [SubNum_f, 1]);
            par_in = fitlm(Profile_All_ext_f, FC_in);
            interval = coefCI(par_in, 0.05);
            interval_value = par_in.Coefficients.Estimate(:) - interval(:,1);
            %[beta p-value] for group
            FC_test_f{i, j} = [par_in.Coefficients.Estimate(4) interval_value(4) par_in.Coefficients.pValue(4)]; 
        else
            FC_test_f{i, j} = nan;
        end
    end    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw heat map of somatosensory regions and thalamus FC
% find soma ROI index and name (onlt include somatosensory function regions)
ROI_soma_index = [];
ROI_soma_name = [];
for i = 1:length(names)
    if  ~isempty(cell2mat(strfind(names(i), 'PostCG'))) || ~isempty(cell2mat(strfind(names(i), 'Parietal Oper'))) || ~isempty(cell2mat(strfind(names(i), 'SPL'))) || ~isempty(cell2mat(strfind(names(i), 'aSMG'))) || ~isempty(cell2mat(strfind(names(i), 'pSMG')))
        ROI_soma_index = [ROI_soma_index, i];
        soma_name_end_u = cell2mat(strfind(names(i), '('));
        soma_name = char(names(i));
        ROI_soma_name = [ROI_soma_name, string(soma_name(15:end))];
    end    
end 

% find only thalamus
ROI_tha_index = [];
ROI_tha_name = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'Thalamus')))
        ROI_tha_index = [ROI_tha_index, i];
        thalamus_name = char(names(i));
        ROI_tha_name = [ROI_tha_name, string(thalamus_name(15:end))];
    end    
end

ROI_tha_name = ["Thalamus R" "Thalamus L"];
ROI_soma_name = ["PostCG R" "PostCG L"	"SPL R" "SPL L" "aSMG R" "aSMG L"	"pSMG R" "pSMG L" "PO R" "PO L"];

Test_ThaSoma_m = sub_extract(ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_m);
Test_ThaSoma_f = sub_extract(ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_f);

%% draw heat map from auditory cortex (ref: Linke, 2018) to salient network/thalamus
ROI_au_index = [];
ROI_au_name = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'HG'))) || ~isempty(cell2mat(strfind(names(i), 'pSTG'))) || ~isempty(cell2mat(strfind(names(i), 'aSTG'))) 
        Au_name = char(names(i));
        Au_name_st = cell2mat(strfind(names(i), '.'));
        Au_name_end = cell2mat(strfind(names(i), '('));
        Au_name_r = cell2mat(strfind(names(i), 'r'));
        Au_name_l = cell2mat(strfind(names(i), 'l'));
        if length(Au_name_st) == 1 
            if Au_name_r(1) < Au_name_end
                ROI_au_index = [ROI_au_index, i];
                Au_name = [Au_name(Au_name_st+1:Au_name_r(1)-1) 'R'];
                ROI_au_name = [ROI_au_name, string(Au_name)];
            elseif Au_name_l(1) < Au_name_end
                ROI_au_index = [ROI_au_index, i];
                Au_name = [Au_name(Au_name_st+1:Au_name_l(3)-1) 'L'];
                ROI_au_name = [ROI_au_name, string(Au_name)];
            else 
                ROI_au_index = [ROI_au_index, i];
                Au_name = [Au_name(Au_name_st+1:Au_name_end-1)];
                ROI_au_name = [ROI_au_name, string(Au_name)];
            end    
        end
    end     
end

Test_ThaAu_m = sub_extract(ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_m);
Test_ThaAu_f = sub_extract(ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_f);

%% draw heat map from visual cortex to salient network/thalamus
ROI_vi_index = [];
ROI_vi_name = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'Visual'))) 
        Vi_name_st = cell2mat(strfind(names(i), '.'));
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

Test_ThaVi_m = sub_extract(ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_m);
Test_ThaVi_f = sub_extract(ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_f);

%% draw heat map from taste cortex (insula) to salient network/thalamus
ROI_ins_index = [];
ROI_ins_name = [];
for i = 1:length(names)
    if ~isempty(cell2mat(strfind(names(i), 'Insular'))) 
        Ins_name_st = cell2mat(strfind(names(i), '.'));
        Ins_name_end = cell2mat(strfind(names(i), '('));
        Ins_name = char(names(i));
        ROI_ins_index = [ROI_ins_index, i];
        ROI_ins_name = [ROI_ins_name, string(Ins_name(Ins_name_st+1:Ins_name_end-1))];     
    end     
end

ROI_ins_name = ["IC R" "IC L"];

Test_ThaIns_m = sub_extract(ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_m);
Test_ThaIns_f = sub_extract(ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_f);

