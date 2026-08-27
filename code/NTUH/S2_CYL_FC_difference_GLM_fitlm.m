%  grouping FC by ROI*ROI into 3D array
clear all
dataDir = ['/Volumes/Milk/CYL/batch_output_0.01_0.08/'];

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

%% read Profile of age, meanFD and gender for the covariates matrix
Profile_All = readcell('/Volumes/Milk/CYL/Profile_sheet/All_profile.xlsx');
% extract age and meanFD only, and gender and FSIQ
Profile_All_ext = [];
for i = 1:SubNum
    Profile_All_ext(i, 1) = str2double(Profile_All{i+1, 2}(1:5));% age
    Profile_All_ext(i, 2) = Profile_All{i+1, 3};% meanFD
    Profile_All_ext(i, 3) = Profile_All{i+1, 4};% sex
    Profile_All_ext(i, 4) = Profile_All{i+1, 5};% Group for GLM test
    Profile_All_ext(i, 5) = Profile_All_ext(i, 4)*Profile_All_ext(i, 3);%sex group interaction
end   

%% assign Z matrix from output into 3D array by group
% NAC group
FC_NAC = [];
% for i = 1:length(files_NAC)
for i = 1:length(files_NAC)
    data_name = string(files_NAC(i));
    disp(data_name);
    ToOpen = append(dataDir, data_name, '/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat');
    load(ToOpen);
    FC_NAC = cat(3, FC_NAC, Z);
end    

% ASD group
FC_ASD = [];
% for i = 1:length(files_ASD)
for i = 1:length(files_ASD)
    data_name = string(files_ASD(i));
    disp(data_name);
    ToOpen = append(dataDir, data_name, '/results/firstlevel/SBC_01/resultsROI_Subject001_Condition001.mat');
    load(ToOpen);
    % Z is the functional connectivity of each functional pair (regionA and
    % regionB), store as a 171*172 z value matrix
    FC_ASD = cat(3, FC_ASD, Z);
end   

% stake up two group to regress out to the same extent
FC_All = cat(3, FC_NAC, FC_ASD);

%% GLM with age, meanFD, gender and Group to test beta 
FC_test = {};
FC_test_sex = {};
FC_test_interact = {};
[dim1, dim2] = size(Z);
for i = 1:dim1
    for j = 1:dim2
        if i ~= j
            FC_in = reshape(FC_All(i,j, :), [SubNum, 1]); 
            par_in = fitlm(Profile_All_ext, FC_in);
            interval = coefCI(par_in, 0.05);
            interval_value = par_in.Coefficients.Estimate(:) - interval(:,1);
            %[beta p-value interval_value] for group
            FC_test{i, j} = [par_in.Coefficients.Estimate(5) interval_value(5) par_in.Coefficients.pValue(5)]; 
            %[beta p-value] for sex
            FC_test_sex{i, j} = [par_in.Coefficients.Estimate(4) interval_value(4) par_in.Coefficients.pValue(4)];
            %[beta p-value] for sex group interaction
            FC_test_interact{i, j} = [par_in.Coefficients.Estimate(6) interval_value(6) par_in.Coefficients.pValue(6)];
        else
            FC_test{i, j} = nan;
            FC_test_sex{i, j} = nan;
            FC_test_interact{i, j} = nan;
        end
    end    
end 

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

Test_ThaSoma = sub_extract(ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test);
Test_sex_ThaSoma = sub_extract(ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_sex);
Test_interact_ThaSoma = sub_extract(ROI_tha_index, ROI_soma_index, ROI_tha_name, ROI_soma_name, FC_test_interact);

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

Test_ThaAu = sub_extract(ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test);
Test_sex_ThaAu = sub_extract(ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_sex);
Test_interact_ThaAu = sub_extract(ROI_tha_index, ROI_au_index, ROI_tha_name, ROI_au_name, FC_test_interact);

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

Test_ThaVi = sub_extract(ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test);
Test_sex_ThaVi = sub_extract(ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_sex);
Test_interact_ThaVi = sub_extract(ROI_tha_index, ROI_vi_index, ROI_tha_name, ROI_vi_name, FC_test_interact);

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

Test_ThaIns = sub_extract(ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test);
Test_sex_ThaIns = sub_extract(ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_sex);
Test_interact_ThaIns = sub_extract(ROI_tha_index, ROI_ins_index, ROI_tha_name, ROI_ins_name, FC_test_interact);
