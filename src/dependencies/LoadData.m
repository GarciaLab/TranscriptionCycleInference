function data_all=LoadData(fileDir,loadPrevious)
% Load data loads data for all operating systems. Directory names
% automatically contain "\" for WINDOWS and "/" for OS,LINUX.

%If loadPrevious setting is specified, choose which results to load.
%Still need to update this section
initialresults_all = [];
if loadPrevious
    initialresults = dir(fullfile(fileDir,'*.mat'));
    initialresults_names = {initialresults.name};

    %Only look at InitialRise data
    initialriseindices = contains(initialresults_names,'InitialRise');
    initialresults = initialresults(initialriseindices);
    initialresults_names = {initialresults_names{initialriseindices}};

    [s,v] = listdlg('PromptString','Select a dataset:','SelectionMode',...
        'multiple','ListString',initialresults_names); %s is the indices of the chosen datasets

    %Load the InitialRise results (just store elongation rate for now, along with Prefix)
    initialresults_all = struct('MCMCresults',{},'Prefix',{}); %Structure containing InitialRise results
    for i = 1:length(s)
        initialresults_temp = load(fullfile(initialresults(s(i)).folder,initialresults_names{s(i)}));
        initialresults_all(i).Prefix = initialresults_temp.Prefix;
        for j = 1:length(initialresults_temp.MCMCresults)
            initialresults_all(i).MCMCresults(j).mean_v = initialresults_temp.MCMCresults(j).mean_v;
            initialresults_all(i).MCMCresults(j).cell_index = initialresults_temp.MCMCresults(j).cell_index;
            initialresults_all(i).MCMCresults(j).ApprovedFits = initialresults_temp.MCMCresults(j).ApprovedFits;

        end
    end

    %Load the data
    fileloc = fullfile('S:','Jonathan','Dropbox','2 Color Elongation','Filtered Particles'); %Loads from 'S:\Jonathan\Dropbox\2 Color Elongation\Filtered Particles' on Windows
    data_all = struct('ElapsedTime',{},'Prefix',{},'TwoColorFilteredParticles',{},...
        'nc12',{},'nc13',{},'nc14',{});
    for i = 1:length(initialresults_all)
        dat = load(fullfile(fileloc,[initialresults_all(i).Prefix,'.mat']));
        temp.ElapsedTime = dat.ElapsedTime;
        temp.Prefix = dat.Prefix;
        temp.nc12 = dat.nc12;
        temp.nc13 = dat.nc13;
        temp.nc14 = dat.nc14;
        temp.TwoColorFilteredParticles = dat.TwoColorFilteredParticles;
        data_all(i) = temp;
    end
else
    %Otherwise, choose which dataset to infer results for.
    files = dir(fullfile(fileDir,'*.mat'));
    names = {files.name};

    [s,v] = listdlg('PromptString','Select a dataset:','SelectionMode',...
        'multiple','ListString',names); %s is the indices of the chosen datasets

    %Load the data
    data_all = struct('data',{});
    for i = 1:length(s)
        filename_dat = fullfile(files(s(i)).folder,names{s(i)});
        dat = load(filename_dat);
        data_all(i).data = dat.data;
    end
end