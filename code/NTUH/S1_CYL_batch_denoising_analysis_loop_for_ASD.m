%% single subject denoising and 1 level analysis for the resting state fMRI in SENSE/VMR project
%% addspm12 first!!!

% create datalist for analysis
clear all
dataDir = ['/Volumes/Milk/CYL/RawData_2_for_ASD']; % find data directory
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
for subNum = 71:length(useable_data)
    data_name = string(useable_data(subNum));
    disp(data_name);
    unsorted_filename = data_name;
    
    % FUNCTIONAL_FILE    
    % assign each file start with 'swuao' functional file in
    search_dir = append('/Volumes/Milk/CYL/RawData_2_for_ASD/', data_name);
    file_name = dir(fullfile(search_dir, 'swau*'));
    file_path = fullfile(search_dir, {file_name.name});
    FUNCTIONAL_FILE = cellstr(file_path');
    
    % STRUCTURAL_FILE
    % take 'wc0c*' as STRUCTURAL_FILE
    search_dir = append('/Volumes/Milk/CYL/RawData_2_for_ASD/', data_name); 
    file_name = dir(fullfile(search_dir, 'wc0c*'));
    file_path = fullfile(search_dir, {file_name.name});
    STRUCTURAL_FILE = cellstr(file_path');

    %% CONN Setup
    cd('/Volumes/Milk/CYL/batch_output_0.01_0.08'); % output file dir
    conn_filename_name = char(append('conn_', data_name, '.mat')); 
    batch.filename=fullfile(pwd, conn_filename_name);
    q_word = append(data_name, ' project already exists, Overwrite?');
    if ~isempty(dir(batch.filename)), 
        Ransw=questdlg(q_word,'warning','Yes','No','No');
        if ~strcmp(Ransw,'Yes')
            continue; 
        end
    end

    batch.Setup.nsubjects=1;
    batch.Setup.structurals{1}= STRUCTURAL_FILE;
    batch.Setup.functionals = cell(1, 1);
    batch.Setup.functionals{1}{1} = FUNCTIONAL_FILE;
    % In our experiment TR == 2
    batch.Setup.RT = 2;
    batch.Setup.conditions.names = {'rest'};
    % initialize resting session
    batch.Setup.conditions.onset{1}{1}{1} = 0;
    batch.Setup.conditions.duration{1}{1}{1} = inf;

    batch.Setup.rois.names={'Default_atlas', 'Default_network', 'Yeo2011'};
    batch.Setup.rois.files{1}=fullfile(fileparts(which('conn')),'rois','atlas.nii');
    batch.Setup.rois.files{2}=fullfile(fileparts(which('conn')),'rois','networks.nii');
    batch.Setup.rois.files{3}={'/Volumes/Milk/CYL/conn/utils/otherrois/Yeo2011.nii'};

    % input the already preprocessed segementation of GM, WM, CSF into mask
    file_name = dir(fullfile(search_dir, '/wc1c*'));
    file_path = fullfile(search_dir, {file_name.name});
    batch.Setup.masks.Grey{1} = cellstr(file_path'); 
    file_name = dir(fullfile(search_dir, '/wc2c*'));
    file_path = fullfile(search_dir, {file_name.name});
    batch.Setup.masks.White{1} = cellstr(file_path'); 
    file_name = dir(fullfile(search_dir, '/wc3c*'));
    file_path = fullfile(search_dir, {file_name.name});
    batch.Setup.masks.CSF{1} = cellstr(file_path'); 

    batch.Setup.covariates.names={'realignment','scrubbing'};
    % assign covariate regressor for resting sessions
    file_name = dir(fullfile(search_dir, '/rp*'));
    file_path = fullfile(search_dir, {file_name.name});
    batch.Setup.covariates.files{1}{1}{1} = cellstr(file_path'); 
    file_name = dir(fullfile(search_dir, '/art_regression_outliers_au*'));
    file_path = fullfile(search_dir, {file_name.name});
    batch.Setup.covariates.files{2}{1}{1} = cellstr(file_path'); 

    % just output ROI to ROI
    batch.Setup.Analysis.type = 1;
    % output BETA map and confound-corrected timeseries
    batch.Setup.outputfiles(1) = 1;
    batch.Setup.outputfiles(2) = 1;
    batch.Setup.isnew=1;
    batch.Setup.done=1;
    
    %% CONN Denoising
    batch.Denoising.filter=[0.01, 0.08];  % frequency filter (band-pass values, in Hz)
    batch.Denoising.confounds.names={'White Matter', 'CSF', 'realignment','scrubbing', 'Effect of rest'}; % Components
    % Dimensions of confound components, scrubbing and effect of rest is
    % variable among subjects, may need excluded for some really off subject
    % check the QC plot in CONN GUI
    batch.Denoising.confounds.dimensions={6 6 6 Inf Inf}; 
    batch.Denoising.done=1;
    
    %% CONN Analysis
    batch.Analysis.analysis_number=1;       % Sequential number identifying each set of independent first-level analyses
    batch.Analysis.measure=1;               % connectivity measure used {1 = 'correlation (bivariate)', 2 = 'correlation (semipartial)', 3 = 'regression (bivariate)', 4 = 'regression (multivariate)';
    batch.Analysis.weight=2;                % within-condition weight used {1 = 'none', 2 = 'hrf', 3 = 'hanning';
    batch.Analysis.sources={};              % (defaults to all ROIs)
    batch.Analysis.done=1;
    
    conn_batch(batch);
    disp(append(data_name, ' done'));
 end