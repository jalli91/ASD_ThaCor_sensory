clear all
dataDir = ['/Volumes/Milk/CYL/RawData_2_for_ASD'];
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
incomplete_data = [];

% for subNum = 1:length(useable_data)
outlier_array = {};
for subNum = 1:length(useable_data)
    data_name = string(useable_data(subNum));
    disp(data_name);
    outlier_array{subNum, 1} = data_name;
    search_dir = append('/Volumes/Milk/CYL/RawData_2_for_ASD/', data_name);
    file_name = dir(fullfile(search_dir, '/art_regression_outliers_au*'));
    ToOpen = append(search_dir, '/', file_name.name);
    load(ToOpen);
    outlier_array{subNum, 2} = sum(sum(R));
    outlier_array{subNum, 3} = 1 - sum(sum(R))/180;
end    