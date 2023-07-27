function data_all=LoadData(fileDir,loadPrevious)
%LoadData loads data for all operating systems. Directory names
%automatically contain "\" for WINDOWS and "/" for OS,LINUX.
% If loadPrevious=true, then the function loads previously processed data
% from dataset files containing 'InitialRise' in their name and
% corresponding dataset files from a fixed location.
% 
% The structures in the InitialRise datasets need to contain the variables
% 'Prefix' and 'MCMCResults'. The variable 'MCMCResults' is a structure
% array containing the fields 'mean_v', 'cell_index' and 'ApprovedFits'.
% The datasets at the fixed location are structures that need to contain
% the fields 'ElapsedTime', 'Prefix', 'TwoColorFilteredParticles', 'nc12',
% 'nc13', and 'nc14'.

%If loadPrevious setting is specified, choose which results to load.
%Still need to update this section
if loadPrevious

    % Specify fixed location
    fileloc = fullfile('S:','Jonathan','Dropbox','2 Color Elongation','Filtered Particles'); %Loads from 'S:\Jonathan\Dropbox\2 Color Elongation\Filtered Particles' on Windows
    % Check if location is already specified
    disp('If you have already specified the dataset location, press {y}, otherwise press any other key.');
    waitforbuttonpress; %User input to press a key
    key = f.CurrentCharacter; %Last pressed key
    if ~strcmp(key,'y')
        fileloc = uigetdir(pwd,'Pick the fixed directory.');
    end

    %Get list of mat-files in fileDir
    initialresults_files = dir(fullfile(fileDir,'*.mat')); % Get structure array of containing file information
    initialresults_names = {initialresults_files.name}; % Extract cell array of file names

    %Only look at InitialRise data (i.e., file names that contain "InitialRise")
    initialriseindices = contains(initialresults_names,'InitialRise'); % Select indices of file names containing "InitialRise"
    initialresults_files = initialresults_files(initialriseindices); % Filter the file list and maintain only those with "InitialRise" in its name
    initialresults_names = initialresults_names(initialriseindices); % Filter the cell array of names

    % Open a dialog box to manually select from the filtered file list
    [s,~] = listdlg('PromptString','Select a dataset:','SelectionMode',...
        'multiple','ListString',initialresults_names); %s is the indices of the chosen datasets

    %Load the InitialRise results (just store elongation rate for now, along with Prefix)
    initialresults_all = struct('MCMCresults',{},'Prefix',{}); %Structure containing InitialRise results
    for i = 1:length(s)
        initialresults_temp = load(fullfile(initialresults_files(s(i)).folder,initialresults_names{s(i)})); %Load previous results
        initialresults_all(i).Prefix = initialresults_temp.Prefix; % Store prefix of dataset i
        for j = 1:length(initialresults_temp.MCMCresults)
            initialresults_all(i).MCMCresults(j).mean_v = initialresults_temp.MCMCresults(j).mean_v; %Store mean elongation rate of nucleus j in dataset i
            initialresults_all(i).MCMCresults(j).cell_index = initialresults_temp.MCMCresults(j).cell_index;
            initialresults_all(i).MCMCresults(j).ApprovedFits = initialresults_temp.MCMCresults(j).ApprovedFits;

        end
    end

    %Load the related datasets "Prefix.mat" from a fixed location fileloc
    data_all = struct('ElapsedTime',{},'Prefix',{},'TwoColorFilteredParticles',{},...
        'nc12',{},'nc13',{},'nc14',{});
    for i = 1:length(initialresults_all)
        dat = load(fullfile(fileloc,[initialresults_all(i).Prefix,'.mat'])); %Load dataset with name [Prefix,'.mat']
        % Temporally store relevant data of dataset i in the structure temp
        temp.ElapsedTime = dat.ElapsedTime;
        temp.Prefix = dat.Prefix;
        temp.nc12 = dat.nc12;
        temp.nc13 = dat.nc13;
        temp.nc14 = dat.nc14;
        temp.TwoColorFilteredParticles = dat.TwoColorFilteredParticles;
        
        % Store temp in data_all(i)
        data_all(i) = temp;
    end
else
    %Otherwise, choose which dataset to infer results for.
    files = dir(fullfile(fileDir,'*.mat'));
    names = {files.name};
    
    % Open a dialog box to manually select from the file list
    [s,~] = listdlg('PromptString','Select a dataset:','SelectionMode',...
        'multiple','ListString',names); %s is the indices of the chosen datasets

    %Load the data
    data_all = struct('data',{});
    for i = 1:length(s)
        filename_dat = fullfile(files(s(i)).folder,names{s(i)});
        dat = load(filename_dat);
        data_all(i).data = dat.data;
    end
end