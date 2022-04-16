function ApproveMCMCResults(varargin)

%Analyzes saved MCMC results of 2 color elongation data. The user has the option of
%approving or rejecting the results of each single nucleus fit. The
%approved/rejected results are saved in the same .mat file as the MCMC results.

%The MCMC results should be saved as a .mat file with 2 stuctures, MCMCplot
%and MCMCresults. Each should be of size 1xN, where N is the number of
%single nucleus fits in the dataset. Additionally, if the file has already
%been curated, there will be a third structure, ApprovedFits. If not, this
%script will create that structure. The values of ApprovedFits correspond
%to as follows:
%   1: approved
%   0: uncurated (can be approved if ApproveAll is used in post-analysis)
%   -1: rejected

%Variable input arguments:
%   'DwellModel': analyze dwell time model results
%   'TerminationModel': analyze termination model results
%   'InitialRise': analyze initial rise analysis results
%   'WholeCycle': analyze whole cycle analysis results
%   'LoadPrevious': load ApprovedFits from previous analysis to save into
%   this one
%   'nc13': analyze nc13
%   'nc14': analyze nc14
%   'RawChains': view raw MCMC chains

%% Input arguments
%By default, user select which datasets to load
load_initialrise = false;
load_wholecycle = false;
load_nc13 = false;
load_nc14 = false;
load_meanrate = false;
load_polyrate = false;
load_previous = false;
RawChains = false; %By default, don't look at raw chains.
span = 0.1; %By default, smooth loading rate using 10% of datapoints.

for i=1:length(varargin)
    if strcmpi(varargin{i},'InitialRise')
        load_initialrise = true;
    end
    if strcmpi(varargin{i},'WholeCycle')
        load_wholecycle = true;
    end
    if strcmpi(varargin{i},'nc13')
        load_nc13 = true;
    end
    if strcmpi(varargin{i},'nc14')
        load_nc14 = true;
    end
    if strcmpi(varargin{i},'MeanRate')
        load_meanrate = true;
    end
    if strcmpi(varargin{i},'PolyRate')
        load_polyrate = true;
    end
    if strcmpi(varargin{i},'RawChains')
        RawChains = true;
    end
    if strcmpi(varargin{i},'SmoothRate')
        span = varargin{i+1};
    end
    if strcmpi(varargin{i},'LoadPrevious')
        load_previous = true;
    end
end
%% Load results
files = dir('S:\Jonathan\Dropbox\2 Color Elongation\MCMC Results\DwellModel\*.mat');
names = {files.name};

%Check for analysis type and nuclear cycle.
s1 = []; %Analysis type
s2 = []; %Nuclear cycle
s3 = 1:length(names);% Rate fitting type (default is to fit all timepoints)
analysistype = [];
ncstring = [];
ratetype = [];
%Check for analysis type
if load_initialrise
    s1 = find(contains(names,'InitialRise'));
    analysistype = 'InitialRise';
elseif load_wholecycle
    s1 = find(contains(names,'Whole'));
    analysistype = 'WholeCycle';
end
%Check for which nuclear cycle
if load_nc13
    s2 = find(contains(names,'nc13'));
    ncstring = 'nc13';
elseif load_nc14
    s2 = find(contains(names,'nc14'));
    ncstring = 'nc14';
end
%Check for rate fitting type
if load_meanrate
    s3 = find(contains(names,'MeanRate'));
    ratetype = 'MeanRate';
elseif load_polyrate
    s3 = find(contains(names,'PolyRate'));
    ratetype = 'PolyRate';
end
%Combine all constraints
s_all = intersect(intersect(s1,s2),s3);
if isempty(s_all)
    s_all = 1:length(names);
end

%Choose which results to load
[s,v] = listdlg('PromptString','Select a dataset:','SelectionMode',...
'single','ListString',names(s_all),'ListSize',[360 300]); %s is the index of the chosen results

%Load data
m = matfile([files(s_all(s)).folder,'\',char(names(s_all(s)))],...
    'Writable',true); %MCMC results .mat file
vars = whos(m); %Unpack the variables
varnames = {vars.name};
N = length(m.MCMCplot); %Number of nuclei fits in this dataset

%Extract the data into temporary variables.
MCMCplot = m.MCMCplot;
MCMCresults = m.MCMCresults;
if RawChains
    chain_name = char(names(s_all(s)));
    chains = load([files(s_all(s)).folder,'\Raw Chains\',chain_name(1:end-4),'_RawChain.mat']);
    MCMCchain = chains.MCMCchain;
end

%If loading previous ApprovedFits, choose file to load
if load_previous
    [file,path] = uigetfile('S:\Jonathan\Dropbox\2 Color Elongation\MCMC Results\DwellModel\*.mat');
    prev = load([path,file]);
    MCMCresults_previous = prev.MCMCresults;
    if length(MCMCresults) ~= length(MCMCresults_previous)
        error('Please select previous analysis of the same dataset');
    end
    for i = 1:length(MCMCresults_previous)
        MCMCresults(i).ApprovedFits = MCMCresults_previous(i).ApprovedFits;
    end
end

%% Approve/reject fits

%User keypress options
InputKey = {'a';'r';',';'.';'j';'x'};
Function = {'Approve this fit';'Reject this fit';'Previous nucleus';...
    'Next nucleus';'Jump to a specific nucleus';'Exit and save'};
t = table(InputKey,Function);

%Figure color map
colormap = {'red',[0.94 0.94 0.94],'green'}; %Approved/Uncurated/Rejected colors

running = true; %Keep running program until stopped
i = 1; %Start with first indexed nucleus for now

close all

if RawChains
    f = figure('Name','Inference Results','Position',[50 100 500 600]);
    c = figure('Name','Parameter distributions','Position',[600 100 600 600]);
else
    f = figure('Name','Inference Results','Position',[200 100 800 600]);
    
end
while running
    
clf(f); %Clear figure
if RawChains
    clf(c);
end

%Extract plotting variables from results
t_plot = MCMCplot(i).t_plot;
MS2_plot = MCMCplot(i).MS2_plot;
PP7_plot = MCMCplot(i).PP7_plot;
%t_interp = MCMCplot(i).t_interp;
%MS2_interp = MCMCplot(i).MS2_interp;
%PP7_interp = MCMCplot(i).PP7_interp;
simMS2 = MCMCplot(i).simMS2;
simPP7 = MCMCplot(i).simPP7;
MeanAP = MCMCplot(i).MeanAP;

%Get fit results
mean_v = MCMCresults(i).mean_v;
sigma_v = MCMCresults(i).sigma_v;
mean_R0 = MCMCresults(i).mean_R(1);
sigma_R0 = MCMCresults(i).sigma_R(1);
mean_dR = MCMCresults(i).mean_R(2:end);
sigma_dR = MCMCresults(i).sigma_R(2:end);
mean_ton = MCMCresults(i).mean_ton;
sigma_ton = MCMCresults(i).sigma_ton;
mean_dwelltime = MCMCresults(i).mean_dwell;
sigma_dwelltime = MCMCresults(i).sigma_dwell;

%Smoothed loading rate
rate_plot = mean_R0 + mean_dR;
rateerror_plot = sigma_R0 + sigma_dR;
ratesmooth_plot = smooth(rate_plot,span);

%Remove loading rates before inferred initiation time
remove_times = find(t_plot(1:end-1) < mean_ton);
rate_plot(remove_times) = nan;
rateerror_plot(remove_times) = nan;
ratesmooth_plot(remove_times) = nan;

%Raw chain results
if RawChains
    dwelltime_chain = MCMCchain(i).dwell_chain; %termination dwell time
    R0_chain = MCMCchain(i).R_chain(:,1); %mean loading rate
    dR_chain = MCMCchain(i).R_chain(:,2:end); %last loading rate fluctuation 
end

%% Plot fit results

subplot(2,1,1)

hold on
plot(t_plot,MS2_plot,'ro','DisplayName','MS2 data');
plot(t_plot,PP7_plot,'go','DisplayName','PP7 data');
%plot(t_interp,MS2_interp,'r-','DisplayName','MS2 interpolation');
%plot(t_interp,PP7_interp,'g-','DisplayName','PP7 interpolation');
plot(t_plot,simMS2,'k--','DisplayName','MS2 fit');
plot(t_plot,simPP7,'k--','DisplayName','PP7 fit');


xlim([t_plot(1), t_plot(end)*1.3]);
ylim([0,  max(PP7_plot) * 1.2]);
xlabel('Time since nuclear cycle start (min)');
ylabel('Fluorescence (AU)');
legend('Location','Northwest')
title({['Single nucleus fit ',num2str(i),' of ',num2str(N),', AP Position ',num2str(MeanAP)]...
    ['Elongation rate = ',num2str(mean_v), ' +/- ',num2str(sigma_v), 'kb/min']});
f.Color = colormap{MCMCresults(i).ApprovedFits+2}; %Set color depending on approval status

%Plot inferred rate
subplot(2,1,2);
hold on
errorbar(t_plot(1:end-1),rate_plot,rateerror_plot,'r.-','CapSize',0,...
    'DisplayName','Inferred loading rate');
plot(t_plot(1:end-1),ratesmooth_plot,'k--','DisplayName','Smoothed loading rate');
line(xlim,mean_R0(1)*[1,1],'Color','black','LineStyle','-','DisplayName',...
    'Inferred mean loading rate');
hold off

xlim([t_plot(1), t_plot(end)*1.3]);
%ylim([0,  max(MS2_plot) * 1.2]);
line(mean_ton*[1,1],ylim,'Color','blue','LineStyle','--',...
    'DisplayName','Inferred time on');
xlabel('Time since nuclear cycle start (min)');
ylabel('Loading rate (AU/min)');
legend('Location','southeast');

%Plot raw chains if desired
if RawChains
    figure(c);
    hold on;

    subplot(3,2,1)
    histogram(dwelltime_chain);
    xlabel('Dwell time (min)');

    subplot(3,2,2)
    plot(dwelltime_chain,'b.');
    ylabel('Dwell time(min)');

    subplot(3,2,3)
    histogram(R0_chain);
    xlabel('Mean loading rate (AU/min)');

    subplot(3,2,4)
    plot(R0_chain,'b.');
    ylabel('Mean loading rate (AU/min)');

    subplot(3,2,5)
    histogram(dR_chain(:,end));
    xlabel('Last loading rate fluctuation (AU/min)');

    subplot(3,2,6)
    plot(dR_chain(:,end),'b.');
    ylabel('Last loading rate fluctuation (AU/min)');
end

% User options (approve/reject, change nucleus)
disp(t); %Display options
exitflag = false; %Loop keypress query until valid exit keypress

while ~exitflag
figure(f);
waitforbuttonpress; %User input to press a key
key = f.CurrentCharacter; %Last pressed key

if strcmp(key,'a')
    MCMCresults(i).ApprovedFits = 1;
    disp('Approved');
elseif strcmp(key,'r')
    MCMCresults(i).ApprovedFits = -1;
    disp('Rejected');
elseif strcmp(key,',')
    if i > 1
        i = i - 1;
        disp('Switching to previous nucleus');
        exitflag = true;
    elseif i == 1
        disp('Already at first nucleus!');
    end
elseif strcmp(key,'.')
    if i < N
        i = i + 1;
        disp('Switching to next nucleus');
        exitflag = true;
    elseif i == N
        disp('Already at last nucleus!');
    end
elseif strcmp(key,'j')
    j = input('Enter in nucleus number to jump to:');
    if j >= 1 && j <= N
        i = j;
        disp(['Switching to nucleus ',num2str(i)]);
        exitflag = true;
    else
        disp(['Error: please enter an integer between 1 and ',num2str(N)]);
    end
elseif strcmp(key,'x')
    disp('Exiting and saving results')
    exitflag = true;
    running = 0;
end
f.Color = colormap{MCMCresults(i).ApprovedFits+2}; %Set color depending on approval status
end

end

%Save approved fits results
m.MCMCresults = MCMCresults;
disp('Results saved.');

close all;

end
