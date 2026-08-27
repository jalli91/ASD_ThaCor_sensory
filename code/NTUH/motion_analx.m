% create datalist for analysis
clear all
dataDir = ['/Volumes/Sea/CYL/ASD_ThaCor/NTUH/RawData_2_for_ASD']; % find data directory
% Get list of all directories and files in the specified folder
file_list = dir(dataDir);

% Initialize an empty cell array to store the matching folder names
useable_data = {};

% Loop through each entry and check if it's a folder and has 4 digits
for i = 1:length(file_list)
    if file_list(i).isdir  % Check if it's a directory
        folder_name = file_list(i).name;
        if length(folder_name) == 4
            useable_data = [useable_data, folder_name];
        end    
    end
end

%  grouping FC by ROI*ROI into 3D array
dataDir = ['/Volumes/Sea/CYL/ASD_ThaCor/NTUH/batch_output_0.01_0.08'];

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

% grouping file names by gender
% Define the file path of gender according indexing file
filePath = '/Volumes/Sea/CYL/ASD_ThaCor/NTUH/Profile_sheet/gender_index.xlsx';
% Read the Excel file into a table
gender_table = readtable(filePath);
gender_table = table2cell(gender_table);

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

% collect mean FD into array - ASD
motion_file_ASD = {};

% for subNum = 1:length(dataname)
for subNum = 1:length(files_ASD)
    data_name = string(files_ASD{subNum}(6:9));
    search_dir = append('/Volumes/Sea/CYL/ASD_ThaCor/NTUH/RawData_2_for_ASD/', data_name);
    data_file = dir(fullfile(search_dir, '/rp*')); 
    data_path = [data_file.folder '/' data_file.name];

    % Load the motion parameters from the file
    motion = load(data_path);
    % Define the conversion factor for rotations (in mm per radian)
    r = 50;
    % Compute differences between consecutive time points
    % This will result in one fewer row than the original data
    diff_motion = diff(motion, 1, 1);
    % Convert rotational differences (columns 4 to 6) from radians to mm
    diff_motion(:,4:6) = diff_motion(:,4:6) * r;
    % Compute FD for each frame as the sum of absolute differences
    FD = sum(abs(diff_motion), 2);
    % Compute the mean FD across all frames (frames 2 to N)
    meanFD = mean(FD);
    % Display the result
    disp(data_name);
    fprintf('Mean Framewise Displacement: %.4f mm\n', meanFD);
    % Add name and meanFD into array
    char_name = char(files_ASD(subNum));
    motion_file_ASD{subNum, 1} = char_name(6:9);
    motion_file_ASD{subNum, 2} = meanFD;

    % obtain scrubbed scan number
    data_file = dir(fullfile(search_dir, '/art_regression_outliers_au*'));
    data_path = [data_file.folder '/' data_file.name];
    load(data_path);
    scrub_scans = sum(sum(R));
    motion_file_ASD{subNum, 3} = scrub_scans;
    motion_file_ASD{subNum, 4} = scrub_scans/180;
end

% collect mean FD into array - ASD
motion_file_NAC = {};

% for subNum = 1:length(dataname)
for subNum = 1:length(files_NAC)
    data_name = string(files_NAC{subNum}(6:9));
    search_dir = append('/Volumes/Sea/CYL/ASD_ThaCor/NTUH/RawData_2_for_ASD/', data_name);
    data_file = dir(fullfile(search_dir, '/rp*')); 
    data_path = [data_file.folder '/' data_file.name];

    % Load the motion parameters from the file
    motion = load(data_path);
    % Define the conversion factor for rotations (in mm per radian)
    r = 50;
    % Compute differences between consecutive time points
    % This will result in one fewer row than the original data
    diff_motion = diff(motion, 1, 1);
    % Convert rotational differences (columns 4 to 6) from radians to mm
    diff_motion(:,4:6) = diff_motion(:,4:6) * r;
    % Compute FD for each frame as the sum of absolute differences
    FD = sum(abs(diff_motion), 2);
    % Compute the mean FD across all frames (frames 2 to N)
    meanFD = mean(FD);
    % Display the result
    disp(data_name);
    fprintf('Mean Framewise Displacement: %.4f mm\n', meanFD);
    % Add name and meanFD into array
    char_name = char(files_NAC(subNum));
    motion_file_NAC{subNum, 1} = char_name(6:9);
    motion_file_NAC{subNum, 2} = meanFD;

    % obtain scrubbed scan number
    data_file = dir(fullfile(search_dir, '/art_regression_outliers_au*'));
    data_path = [data_file.folder '/' data_file.name];
    load(data_path);
    scrub_scans = sum(sum(R));
    motion_file_NAC{subNum, 3} = scrub_scans;
    motion_file_NAC{subNum, 4} = scrub_scans/180;
end

FD_all = [motion_file_ASD{:,2} motion_file_NAC{:,2}];
disp('ALL mean/std');
mean_FD_all = mean(FD_all);
std_FD_all = std(FD_all);
disp([mean_FD_all std_FD_all]);

disp('NAC mean/std');
disp([mean([motion_file_NAC{:,2}]) std([motion_file_NAC{:,2}])]);
disp('ASD mean/std');
disp([mean([motion_file_ASD{:,2}]) std([motion_file_ASD{:,2}])]);
[h_FD, p_FD] = ttest2([motion_file_ASD{:,2}], [motion_file_NAC{:,2}]);
[h_sc, p_sc] = ttest2([motion_file_ASD{:,4}], [motion_file_NAC{:,4}]);

% mark "off" if subject FD > 2 folds of total standard deviation
for i = 1:length(files_ASD)
    if motion_file_ASD{i, 2} - mean_FD_all > 2*std_FD_all
        motion_file_ASD{i, 5} = "off";
    end    
end

for i = 1:length(files_NAC)
    if motion_file_NAC{i, 2} - mean_FD_all > 2*std_FD_all
        motion_file_NAC{i, 5} = "off";
    end    
end

% gender subgrouping comparsion
motion_file_ASD_m = {};
j = 1;
for i = 1:length(files_ASD)
    if ismember(files_ASD{1, i}, files_ASD_m)
        for n = 1:4
            motion_file_ASD_m{j, n} = motion_file_ASD{i, n};
        end
        j = j + 1;
    end    
end   

motion_file_ASD_f = {};
j = 1;
for i = 1:length(files_ASD)
    if ismember(files_ASD{1, i}, files_ASD_f)
        for n = 1:4
            motion_file_ASD_f{j, n} = motion_file_ASD{i, n};
        end
        j = j + 1;
    end    
end 

motion_file_NAC_m = {};
j = 1;
for i = 1:length(files_NAC)
    if ismember(files_NAC{1, i}, files_NAC_m)
        for n = 1:4
            motion_file_NAC_m{j, n} = motion_file_NAC{i, n};
        end
        j = j + 1;
    end    
end 

motion_file_NAC_f = {};
j = 1;
for i = 1:length(files_NAC)
    if ismember(files_NAC{1, i}, files_NAC_f)
        for n = 1:4
            motion_file_NAC_f{j, n} = motion_file_NAC{i, n};
        end
        j = j + 1;
    end    
end 

[h_FD_m, p_FD_m] = ttest2([motion_file_ASD_m{:,2}], [motion_file_NAC_m{:,2}]);
[h_sc_m, p_sc_m] = ttest2([motion_file_ASD_m{:,4}], [motion_file_NAC_m{:,4}]);

[h_FD_f, p_FD_f] = ttest2([motion_file_ASD_f{:,2}], [motion_file_NAC_f{:,2}]);
[h_sc_f, p_sc_f] = ttest2([motion_file_ASD_f{:,4}], [motion_file_NAC_f{:,4}]);
