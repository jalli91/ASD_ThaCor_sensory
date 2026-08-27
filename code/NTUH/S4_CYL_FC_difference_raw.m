%  grouping FC by ROI*ROI into 3D array
clear all
dataDir = ['/Volumes/Sea/CYL/ASD_ThaCor/NTUH/batch_output_0.01_0.08/'];

%% grouping the result folder in the directory
file_list = dir(dataDir);
files_NAC = {};
files_ASD = {};
for i = 1:length(file_list)
    file_name = file_list(i).name;
    if startsWith(file_name, 'conn_1') & ~endsWith(file_name, '.mat')
        files_NAC{end+1} = file_name;
    elseif startsWith(file_name, 'conn_2') & ~endsWith(file_name, '.mat')
        files_ASD{end+1} = file_name;
    end
end
SubNum = length(files_NAC)+length(files_ASD);

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
    FC_ASD = cat(3, FC_ASD, Z);
end   
FC_All = cat(3, FC_NAC, FC_ASD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find ROI index and name
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract raw FC of Thalamus-sensory
%% ThaAu
sub_ThaAu_NAC = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_au_name)
        sub_ThaAu_NAC{1, counter} = ROI_tha_name(i) + '*' + ROI_au_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_NAC
for subNum = 1:size(FC_NAC, 3)
    char_name = char(files_NAC(subNum));
    sub_ThaAu_NAC{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_au_index)
            sub_ThaAu_NAC{subNum+1, counter} = FC_NAC(ROI_tha_index(i), ROI_au_index(j), subNum);
            counter = counter + 1;
        end
    end
end 

sub_ThaAu_ASD = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_au_name)
        sub_ThaAu_ASD{1, counter} = ROI_tha_name(i) + '*' + ROI_au_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_ASD
for subNum = 1:size(FC_ASD, 3)
    char_name = char(files_ASD(subNum));
    sub_ThaAu_ASD{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_au_index)
            sub_ThaAu_ASD{subNum+1, counter} = FC_ASD(ROI_tha_index(i), ROI_au_index(j), subNum);
            counter = counter + 1;
        end
    end
end 

%% ThaSoma
sub_ThaSoma_NAC = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_soma_name)
        sub_ThaSoma_NAC{1, counter} = ROI_tha_name(i) + '*' + ROI_soma_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_NAC
for subNum = 1:size(FC_NAC, 3)
    char_name = char(files_NAC(subNum));
    sub_ThaSoma_NAC{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_soma_index)
            sub_ThaSoma_NAC{subNum+1, counter} = FC_NAC(ROI_tha_index(i), ROI_soma_index(j), subNum);
            counter = counter + 1;
        end
    end
end 

sub_ThaSoma_ASD = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_soma_name)
        sub_ThaSoma_ASD{1, counter} = ROI_tha_name(i) + '*' + ROI_soma_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_ASD
for subNum = 1:size(FC_ASD, 3)
    char_name = char(files_ASD(subNum));
    sub_ThaSoma_ASD{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_soma_index)
            sub_ThaSoma_ASD{subNum+1, counter} = FC_ASD(ROI_tha_index(i), ROI_soma_index(j), subNum);
            counter = counter + 1;
        end
    end
end

%% ThaVi
sub_ThaVi_NAC = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_vi_name)
        sub_ThaVi_NAC{1, counter} = ROI_tha_name(i) + '*' + ROI_vi_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_NAC
for subNum = 1:size(FC_NAC, 3)
    char_name = char(files_NAC(subNum));
    sub_ThaVi_NAC{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_vi_index)
            sub_ThaVi_NAC{subNum+1, counter} = FC_NAC(ROI_tha_index(i), ROI_vi_index(j), subNum);
            counter = counter + 1;
        end
    end
end 

sub_ThaVi_ASD = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_vi_name)
        sub_ThaVi_ASD{1, counter} = ROI_tha_name(i) + '*' + ROI_vi_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_ASD
for subNum = 1:size(FC_ASD, 3)
    char_name = char(files_ASD(subNum));
    sub_ThaVi_ASD{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_vi_index)
            sub_ThaVi_ASD{subNum+1, counter} = FC_ASD(ROI_tha_index(i), ROI_vi_index(j), subNum);
            counter = counter + 1;
        end
    end
end

%% ThaIns
sub_ThaIns_NAC = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_ins_name)
        sub_ThaIns_NAC{1, counter} = ROI_tha_name(i) + '*' + ROI_ins_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_NAC
for subNum = 1:size(FC_NAC, 3)
    char_name = char(files_NAC(subNum));
    sub_ThaIns_NAC{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_ins_index)
            sub_ThaIns_NAC{subNum+1, counter} = FC_NAC(ROI_tha_index(i), ROI_ins_index(j), subNum);
            counter = counter + 1;
        end
    end
end 

sub_ThaIns_ASD = {};
counter = 2;
for i = 1:length(ROI_tha_name) % row name as ROIa*ROIb
    for j = 1:length(ROI_ins_name)
        sub_ThaIns_ASD{1, counter} = ROI_tha_name(i) + '*' + ROI_ins_name(j);
        counter = counter + 1;
    end    
end
% raw FC from FC_ASD
for subNum = 1:size(FC_ASD, 3)
    char_name = char(files_ASD(subNum));
    sub_ThaIns_ASD{subNum+1, 1} = char_name(6:9); % column name: subject 
    counter = 2;
    for i = 1:length(ROI_tha_index)
        for j = 1:length(ROI_ins_index)
            sub_ThaIns_ASD{subNum+1, counter} = FC_ASD(ROI_tha_index(i), ROI_ins_index(j), subNum);
            counter = counter + 1;
        end
    end
end