function ApproveMCMCResults(varargin)

%Analyzes saved MCMC results of 2 color elongation data of a particular
%construct (specified in the library). The user has the option of approving
%or rejecting the results of each single nucleus fit. The approved/rejected
%results are saved in the same .mat file as the MCMC results.
%
%The MCMC results should be saved as a .mat file with 2 stuctures, MCMCplot
%and MCMCresults. Each should be of size 1xN, where N is the number of
%single nucleus fits in the dataset. Additionally, if the file has already
%been curated, there will be a third structure, ApprovedFits. If not, this
%script will create that structure. The values of ApprovedFits correspond
%to as follows:
%   1: approved
%   0: uncurated (can be approved if ApproveAll is used in post-analysis)
%   -1: rejected
%
% Variable input arguments:
%   'fileDir': specify file directory. Otherwise a dialog box will open to
%              let you choose the directory.
%   'RawChains': view raw MCMC chains (show default parameter selection)
%   'RawChainsCheckbox': view raw MCMC chains and select the parameters via
%                        a dialog box
%   'InitiationFluctuations': option to simultaneously view the initiation
%                             rate fluctuations 
%   'SmoothRate': percentage of datapoints used for a smoothed plot of the
%                           initiation rate fluctuations 
%   'InitialRise': option to reduce import selection to files containing
%                  "InitialRise" in their name
%
%   'MeanRate': option to reduce import selection to files containing
%               "MeanRate" in their name (inhibits 'PolyRate' option)
%   'PolyRate': option to reduce import selection to files containing
%               "PolyRate" in their name
%   'WholeCycle': option to reduce import selection to files containing
%                 "Whole" in their name
%   'LoadPrevious': import the results of a previous data curation
%                   (ApprovedFits) and overwrite the respective structures
%                   in the currently loaded dataset 
%   'nc13': option to reduce import selection to files containing
%           "nc13" in their name (inhibits 'nc14' option)
%   'nc14': option to reduce import selection to files containing
%           "nc14" in their name

%% Input arguments
%By default, user select which datasets to load
load_initialrise = false;
load_wholecycle = false;
load_nc13 = false;
load_nc14 = false;
load_meanrate = false;
load_polyrate = false;
load_previous = false;

%By default, don't look at raw chains.
RawChains = false;
RawChainsCheckbox = false;


%By default, don't look at initiation fluctuations.
InitiationFluctuations = false;
span = 0.1; %By default, smooth loading rate using 10% of datapoints.
fileDir = '';

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
    if strcmpi(varargin{i},'RawChainsCheckbox')
        RawChainsCheckbox = true;
    end
    if strcmpi(varargin{i},'InitiationFluctuations')
        InitiationFluctuations = true;
    end
    if strcmpi(varargin{i},'SmoothRate')
        span = varargin{i+1};
    end
    if strcmpi(varargin{i},'LoadPrevious')
        load_previous = true;
    end
    if strcmpi(varargin{i},'fileDir')
        fileDir = varargin{i+1};
    end
end

%% Load results
if isempty(fileDir)
    msg_box = msgbox('Choose the directory, that contains the MCMC results.', 'Directory Selection'); % Display a message box with instructions
    uiwait(msg_box); % Halt execution (make the figure modal) until the OK button is clicked on the message box
    fileDir = uigetdir(pwd); % Open dialog box to choose the directory
end
files = dir(fileDir);
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
[s,~] = listdlg('PromptString','Select a dataset:','SelectionMode',...
    'single','ListString',names(s_all),'ListSize',[360 300]); %s is the index of the chosen results

%Load data: Create a matfile object of the chosen dataset (MCMC results). This object
%allows for accessing and changing parts of a .mat file stored on disk
%without loading the whole file into memory.
m = matfile(fullfile(files(s_all(s)).folder,char(names(s_all(s))),filesep),...
    'Writable',true);
% vars = whos(m); %Unpack the variables of m
% varnames = {vars.name};
N = length(m.MCMCplot); %Number of nuclei fits in this dataset

%Extract the data into temporary variables.
MCMCplot = m.MCMCplot;
MCMCresults = m.MCMCresults;

% In case of RawChain option, look for *_RawChain.mat file in fileDir and
% fileDir/Raw Chains/. Otherwise open dialog box to choose the file.
if or(or(RawChains,RawChainsCheckbox),InitiationFluctuations)
    chain_file = fullfile(files(s_all(s)).folder,'Raw Chains',[names{s_all(s)}(1:end-4),'_RawChain.mat']);
    if exist(chain_file, 'file') == 2
        %disp('File exists')
        chains = load(chain_file);
    else
        chain_file = fullfile(files(s_all(s)).folder,[names{s_all(s)}(1:end-4),'_RawChain.mat']);
        if exist(chain_file, 'file') == 2
            %disp('File exists')
            chains = load(chain_file);
        else
            %disp('File does not exist')
            msg_box = msgbox(['Choose the mat-file containing the raw chains assciated with ',names{s_all(s)}(1:end-4),'.'], 'File Selection'); % Display a message box with instructions
            uiwait(msg_box); % Halt execution (make the figure modal) until the OK button is clicked on the message box
            [chain_file,chain_path] = uigetfile(fullfile(files(s_all(s)).folder, '*_RawChain.mat'));
            if contains(chain_file,names{s_all(s)}(1:end-4))
                chains = load([chain_path,chain_file]);
            else
                error('You picked the wrong file for the RawChains or InitiationFluctuations option.');
            end
        end
    end
    % Extract the chains into temporary variable.
    MCMCchain = chains.MCMCchain;
end

%If loading previous ApprovedFits, choose file to load
if load_previous
    [file,path] = uigetfile(pwd,'Choose the dataset containing information on previous approval.');
    prev = load([path,file]);
    MCMCresults_previous = prev.MCMCresults;
    if length(MCMCresults) ~= length(MCMCresults_previous)
        error('Please select previous analysis of the same dataset.');
    end
    for i = 1:length(MCMCresults_previous)
        MCMCresults(i).ApprovedFits = MCMCresults_previous(i).ApprovedFits;
    end
end

%% Load construct details

% Extract construct name from dataset (assuming that all cells in the set
% are of the same kind)
construct = MCMCresults(1).FittedConstruct;

%Query to the construct library regarding the construct details:
%'segments' and 'velocities' contain details about the elongation rates on
%different segments of the construct and the length of the construct
[ElongationSegments,~,~] = library(construct);

%Extract distinguishable names of velocity parameters
velocity_names = unique(ElongationSegments.velocities);

% Generate field names of velocity parameters
velocity_fields_MCMCresults = [cellfun(@(x) strcat('mean_',x), velocity_names, 'UniformOutput', false),...
    cellfun(@(x) strcat('CI_',x), velocity_names, 'UniformOutput', false)];

if RawChainsCheckbox
    % Select chains via dialog box
    [disp_default_chains,defaultText,disp_velocity_chains] = CheckboxDialog(velocity_names);
    
    % Get number of chains to view
    size_rawchainplot = sum(disp_default_chains) + sum(disp_velocity_chains); %number of raw chains displayed (number of subplot rows)
    
    % Obtain names of parameters to view
    defaultParams = defaultText(2,:); defaultParams = defaultParams(disp_default_chains);
    defaultOptions = defaultText(1,:); defaultOptions = defaultOptions(disp_default_chains);
    defaultOption_Units = defaultText(3,:); defaultOption_Units = defaultOption_Units(disp_default_chains);
    velocityOptions = velocity_names(disp_velocity_chains);

    % Generate field names of chains to view
    default_fields_MCMCchains = cellfun(@(x) strcat(x,'_chain'), defaultParams, 'UniformOutput', false);
    velocity_fields_MCMCchains = cellfun(@(x) strcat(x,'_chain'), velocity_names, 'UniformOutput', false);
    velocity_fields_MCMCchains = velocity_fields_MCMCchains(disp_velocity_chains);
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
cellNum = 1; %Start with first indexed nucleus for now

close all

if or(or(RawChains,RawChainsCheckbox),InitiationFluctuations)
    % Split the display with the two plots
    % Get the size of the screen
    screenSize = get(groot, 'Screensize');

    % Set up the figure properties
    left = screenSize(1);
    bottom = screenSize(2);
    width = screenSize(3) / 2;  % Split the width by half
    height = screenSize(4);

    % Define the figures
    f = figure('Name','Inference Results','Position',[left bottom width*0.9 height*0.9]);

    if or(RawChains,RawChainsCheckbox)
    c = figure('Name','Parameter distributions','Position',[left+1.1*width bottom width*0.9 height*0.9]);
    end
    % f = figure('Name','Inference Results','Position',[50 100 500 600]);
    % c = figure('Name','Parameter distributions','Position',[600 100 600 600]);
else
    f = figure('Name','Inference Results','Position',[200 100 800 600]);

end
while running

    clf(f); %Clear figure
    if or(RawChains,RawChainsCheckbox)
        clf(c);
    end

    %Extract plotting variables from results
    t_plot = MCMCplot(cellNum).t_plot;
    MS2_plot = MCMCplot(cellNum).MS2_plot;
    PP7_plot = MCMCplot(cellNum).PP7_plot;
    %t_interp = MCMCplot(cellNum).t_interp;
    %MS2_interp = MCMCplot(cellNum).MS2_interp;
    %PP7_interp = MCMCplot(cellNum).PP7_interp;
    simMS2 = MCMCplot(cellNum).simMS2;
    simPP7 = MCMCplot(cellNum).simPP7;
    %MeanAP = MCMCplot(cellNum).MeanAP;

    %Get fit results of elongation rates, R, dR, ton and tau

    N_idx=length(velocity_fields_MCMCresults)/2;
    mean_velocity = zeros(1,N_idx); %Generate vector of mean parameters
    for idx=1:N_idx
        mean_velocity(idx) = MCMCresults(cellNum).(velocity_fields_MCMCresults{idx});   %Extract and save means of the velocities
        CI_velocities(:,idx) = MCMCresults(cellNum).(velocity_fields_MCMCresults{N_idx+idx}); %Extract and save credible interval of the velocities
    end

    mean_R = MCMCresults(cellNum).mean_R;
    CI_R = MCMCresults(cellNum).CI_R;
    mean_dR_end = MCMCresults(cellNum).mean_dR(end);
    CI_dR_end = MCMCresults(cellNum).CI_dR(:,end);
    mean_ton = MCMCresults(cellNum).mean_ton;
    CI_ton = MCMCresults(cellNum).CI_ton;
    mean_dwelltime = MCMCresults(cellNum).mean_tau;
    CI_dwelltime = MCMCresults(cellNum).CI_tau;

    %Raw chain results
    if RawChains
        dwelltime_chain = MCMCchain(cellNum).tau_chain; %termination dwell time
        R_chain = MCMCchain(cellNum).R_chain; %mean loading rate
        dR_chain_end = MCMCchain(cellNum).dR_chain(:,end); %last loading rate fluctuation
        ton_chain = MCMCchain(cellNum).ton_chain;
    end

    if RawChainsCheckbox

        default_chains = cell(4,sum(disp_default_chains));
        for i = 1:sum(disp_default_chains)
            default_chains{1,i} = MCMCchain(cellNum).(default_fields_MCMCchains{i});
            default_chains{2,i} = defaultOptions{i};
            default_chains{3,i} = MCMCresults(cellNum).(['mean_',defaultParams{i}])(:,end);
            default_chains{4,i} = MCMCresults(cellNum).(['CI_',defaultParams{i}])(:,end);
        end

        velocity_chains = cell(4,sum(disp_velocity_chains));
        v_chains_idx = find(disp_velocity_chains);
        for i = 1:sum(disp_velocity_chains)
            velocity_chains{1,i} = MCMCchain(cellNum).(velocity_fields_MCMCchains{v_chains_idx(i)});
            velocity_chains{2,i} = velocity_names{v_chains_idx(i)};
            velocity_chains{3,i} = mean_velocity(v_chains_idx(i));
            velocity_chains{4,i} = CI_velocities(:,v_chains_idx(i));
        end
    end

    % Initiation fluctuations
    if InitiationFluctuations
        R_chain = MCMCchain(cellNum).R_chain; %mean loading rate
        noisyR_chain = MCMCchain(cellNum).R_chain + MCMCchain(cellNum).dR_chain;
        BayesianCoverage = 0.95; % Set Bayesian coverage of credible intervals
        lower_quantile = (1 - BayesianCoverage)/2; upper_quantile = 1 - lower_quantile; % Set quantiles
        noisyR_CIs = quantile(noisyR_chain,[lower_quantile,upper_quantile]);

        %Smoothed loading rate
        rate_plot = mean(noisyR_chain);
        rateerror_plot_lower = rate_plot - noisyR_CIs(1,:);
        rateerror_plot_upper = noisyR_CIs(2,:) - rate_plot;
        ratesmooth_plot = smooth(rate_plot,span);

        %Remove loading rates before inferred initiation time
        remove_times = find(t_plot(1:end-1) < mean_ton);
        rate_plot(remove_times) = NaN;
        rateerror_plot(remove_times) = NaN;
        ratesmooth_plot(remove_times) = NaN;
    end

    %% Plot fit results

    figure(f);
    if InitiationFluctuations
        subplot(2,1,1)
    end
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
    legend('Location','Northwest');
    % Generate title with multiple lines (displaying a list of all elongation rates fit falues)
    title_array = cell(1,N_idx+2);
    title_array{1} = ['Single nucleus fit ',num2str(cellNum),' of ',num2str(N)];
    title_array{2} = '-------------------------------------------------------------------';
    for idx=1:N_idx
        title_array{idx+2} = ['Elongation rate ',velocity_names{idx},' = ',num2str(mean_velocity(idx)),' kb/min', ' (95% CI: [',num2str(CI_velocities(1,idx)),',',num2str(CI_velocities(2,idx)),'])'];
    end
    title(title_array);
    f.Color = colormap{MCMCresults(cellNum).ApprovedFits+2}; %Set color depending on approval status

    if InitiationFluctuations
        %Plot inferred rate
        subplot(2,1,2);
        hold on
        disp(size(rate_plot));
        disp(size(t_plot));
        errorbar(t_plot,rate_plot,rateerror_plot_lower,rateerror_plot_upper,'r.-','CapSize',0,...
            'DisplayName','Inferred loading rate');
        plot(t_plot,ratesmooth_plot,'k--','DisplayName','Smoothed loading rate');
        line(xlim,mean_R*[1,1],'Color','black','LineStyle','-','DisplayName',...
            'Inferred mean loading rate');
        hold off

        xlim([t_plot(1), t_plot(end)*1.3]);
        ylim([min(rate_plot - rateerror_plot_lower)*(1-sign(min(rate_plot - rateerror_plot_lower))*0.1),  max(rateerror_plot_upper +rate_plot) * 1.2]);
        line(mean_ton*[1,1],ylim,'Color','blue','LineStyle','--',...
            'DisplayName','Inferred time on');
        xlabel('Time since nuclear cycle start (min)');
        ylabel('Loading rate (AU/min)');
        legend('Location','southeast');
    end


    %Plot raw chains if desired
    if RawChains
        figure(c);
        xlim('auto');
        hold on;

        subplot(4,2,1)
        histogram(dwelltime_chain);
        line(mean_dwelltime*[1,1],ylim,'Color','green','LineStyle','--');
        line(CI_dwelltime(1)*[1,1],ylim,'Color','black','LineStyle','--');
        line(CI_dwelltime(2)*[1,1],ylim,'Color','black','LineStyle','--');
        xlabel('Dwell time (min)');
        title(['Mean: ',num2str(mean_dwelltime),' (95% CI: [',num2str(CI_dwelltime(1)),', ',num2str(CI_dwelltime(2)),'])']);

        subplot(4,2,2)
        plot(dwelltime_chain,'b.');
        ylabel('Dwell time(min)');

        subplot(4,2,3)
        histogram(R_chain);
        line(mean_R*[1,1],ylim,'Color','green','LineStyle','--');
        line(CI_R(1)*[1,1],ylim,'Color','black','LineStyle','--');
        line(CI_R(2)*[1,1],ylim,'Color','black','LineStyle','--');
        xlabel('Mean loading rate (AU/min)');
        title(['Mean: ',num2str(mean_R),' (95% CI: [',num2str(CI_R(1)),', ',num2str(CI_R(2)),'])']);

        subplot(4,2,4)
        plot(R_chain,'b.');
        ylabel('Mean loading rate (AU/min)');

        subplot(4,2,5)
        histogram(dR_chain_end);
        line(mean_dR_end*[1,1],ylim,'Color','green','LineStyle','--');
        line(CI_dR_end(1)*[1,1],ylim,'Color','black','LineStyle','--');
        line(CI_dR_end(2)*[1,1],ylim,'Color','black','LineStyle','--');
        xlabel('Last loading rate fluctuation (AU/min)');
        title(['Mean: ',num2str(mean_dR_end),' (95% CI: [',num2str(CI_dR_end(1)),', ',num2str(CI_dR_end(2)),'])']);

        subplot(4,2,6)
        plot(dR_chain_end,'b.');
        ylabel('Last loading rate fluctuation (AU/min)');

        subplot(4,2,7)
        histogram(ton_chain);
        xlabel('on-time (min)');
        line(mean_ton*[1,1],ylim,'Color','green','LineStyle','--');
        line(CI_ton(1)*[1,1],ylim,'Color','black','LineStyle','--');
        line(CI_ton(2)*[1,1],ylim,'Color','black','LineStyle','--');
        title(['Mean: ',num2str(mean_ton),' (95% CI: [',num2str(CI_ton(1)),', ',num2str(CI_ton(2)),'])']);

        subplot(4,2,8)
        plot(ton_chain,'b.');
        ylabel('on-time (min)');
    end

    if RawChainsCheckbox
        figure(c);
        xlim('auto');
        hold on;

        for num_subplot = 1: sum(disp_default_chains)   
        subplot(size_rawchainplot,2,2*num_subplot-1)
        histogram(default_chains{1,num_subplot});
        line(default_chains{3,num_subplot}*[1,1],ylim,'Color','green','LineStyle','--');
        line(default_chains{4,num_subplot}(1)*[1,1],ylim,'Color','black','LineStyle','--');
        line(default_chains{4,num_subplot}(2)*[1,1],ylim,'Color','black','LineStyle','--');
        xlabel([defaultOptions{num_subplot},defaultOption_Units{num_subplot}]);
        title(['Mean: ',num2str(default_chains{3,num_subplot}),' (95% CI: [',num2str(default_chains{4,num_subplot}(1)),', ',num2str(default_chains{4,num_subplot}(2)),'])']);

        subplot(size_rawchainplot,2,2*num_subplot)
        plot(default_chains{1,num_subplot},'b.');
        ylabel([defaultOptions{num_subplot},defaultOption_Units{num_subplot}]);
        end
        for num_subplot = 1:sum(disp_velocity_chains)
        subplot(size_rawchainplot,2,2*num_subplot - 1 + 2*sum(disp_default_chains))
        histogram(velocity_chains{1,num_subplot});
        line(velocity_chains{3,num_subplot}*[1,1],ylim,'Color','green','LineStyle','--');
        line(velocity_chains{4,num_subplot}(1)*[1,1],ylim,'Color','black','LineStyle','--');
        line(velocity_chains{4,num_subplot}(2)*[1,1],ylim,'Color','black','LineStyle','--');
        xlabel(['Elongation rate ',velocity_chains{2,num_subplot},' (kb/min)']);
        title(['Mean: ',num2str(velocity_chains{3,num_subplot}),' (95% CI: [',num2str(velocity_chains{4,num_subplot}(1)),', ',num2str(velocity_chains{4,num_subplot}(2)),'])']);

        subplot(size_rawchainplot,2,2*num_subplot + 2*sum(disp_default_chains))
        plot(velocity_chains{1,num_subplot},'b.');
        ylabel(['Elongation rate ',velocity_chains{2,num_subplot},' (kb/min)']);
        end
    end


    % User options (approve/reject, change nucleus)
    disp(t); %Display options
    exitflag = false; %Loop keypress query until valid exit keypress

    while ~exitflag
        figure(f);
        waitforbuttonpress; %User input to press a key
        key = f.CurrentCharacter; %Last pressed key

        if strcmp(key,'a')
            MCMCresults(cellNum).ApprovedFits = 1;
            disp('Approved');
        elseif strcmp(key,'r')
            MCMCresults(cellNum).ApprovedFits = -1;
            disp('Rejected');
        elseif strcmp(key,',')
            if cellNum > 1
                cellNum = cellNum - 1;
                disp('Switching to previous nucleus');
                exitflag = true;
            elseif cellNum == 1
                disp('Already at first nucleus!');
            end
        elseif strcmp(key,'.')
            if cellNum < N
                cellNum = cellNum + 1;
                disp('Switching to next nucleus');
                exitflag = true;
            elseif cellNum == N
                disp('Already at last nucleus!');
            end
        elseif strcmp(key,'j')
            j = input('Enter in nucleus number to jump to:');
            if j >= 1 && j <= N
                cellNum = j;
                disp(['Switching to nucleus ',num2str(cellNum)]);
                exitflag = true;
            else
                disp(['Error: please enter an integer between 1 and ',num2str(N)]);
            end
        elseif strcmp(key,'x')
            disp('Exiting and saving results')
            exitflag = true;
            running = 0;
        end
        f.Color = colormap{MCMCresults(cellNum).ApprovedFits+2}; %Set color depending on approval status
    end

end

%Save approved fits results
m.MCMCresults = MCMCresults;
disp('Results saved.');

close all;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions

function [leftValues,leftOptions, rightValues] = CheckboxDialog(velocity_names)
    % Number of checkboxes in the right column
    numRightCheckboxes = length(velocity_names);

    % Create a list of permanent options
    leftOptions = {'Termination dwell time', 'Mean initiation rate', 'Last initiation rate fluctuation', 'on-time';...
                    'tau','R','dR','ton';...
                    ' (min)',' (AU/min)',' (AU/min)', ' (min)'};

    % Create a figure and panel within it
    d = dialog('Units', 'normalized', 'Position',[0.2 0.2 0.15 0.15],'Name','Select 3 or 4 Parameters');

    % Create checkboxes for the left column
    leftValues = false(4,1);
    for i = 1:4
        uicontrol(d,'Style','checkbox','String',leftOptions{1,i},...
                  'Units', 'normalized', 'Position',[0.1 0.8-(i-1)*0.2 0.4 0.15],...
                  'HandleVisibility','off', 'Callback', {@leftCheckboxCallback, i});
    end

    % Create checkboxes for the right column
    rightValues = false(numRightCheckboxes, 1);
    for i = 1:numRightCheckboxes
        uicontrol(d,'Style','checkbox','String',velocity_names{i},...
                  'Units', 'normalized', 'Position',[0.6 0.8-(i-1)*0.2 0.4 0.15],...
                  'HandleVisibility','off', 'Callback', {@rightCheckboxCallback, i});
    end

    % Create OK pushbutton
    btn = uicontrol('Parent',d,...
               'Units', 'normalized', 'Position',[0.3 0.05 0.4 0.2],...
               'String','OK',...
               'Callback',@buttonCallback);
    uiwait(d); % Halt the execution until the dialog is deleted

    function leftCheckboxCallback(hObject, eventdata, idx)
        % Update the value of the checkbox
        leftValues(idx) = hObject.Value;
    end

    function rightCheckboxCallback(hObject, eventdata, idx)
        % Update the value of the checkbox
        rightValues(idx) = hObject.Value;
    end

    function buttonCallback(hObject, eventdata)
        % Delete the dialog
        delete(d);
    end
end