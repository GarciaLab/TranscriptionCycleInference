function TranscriptionCycleMCMC(varargin)

%Fits transcription cycle model from Liu et al (2020) using Markov Chain
%Monte Carlo, implemented in the MCMCstat package. The data must be
%generated using a dual color MS2/PP7 reported as described in Liu et al.

%The code can load multiple datasets, and possesses user-specified options
%as variable arguments. The details of these inputs are described in the
%readme and in the manuscript, and are initialized with default values if
%not specified by the user.


%Variable arguments:
%   'n_burn', n_burn: number of burn-in steps in MCMC algorithm (default = 10000)
%   'n_steps', n_steps: number of steps in MCMC algorithm including burn-in (default =
%   20000)
%   'ratePriorWidth', ratePriorWidth: standard deviation of Gaussian prior
%   for rate fluctuation term dR(t) (default = 50, see Liu et al, Section S3.1) 
%   't_start', t_start: time value to start fit at
%   't_end', t_end: time value to end fit at
%   'loadPrevious', true/false: option to load previous inference results
%   to retain inferred elongation rate for hierarchical fit (see Liu et al,
%   Section S3.2
%   'construct', construct: option to specify a custom reporter gene that
%   must be defined in the subfunction GetFluorFromPolPos (default uses the
%   P2P-MS2-lacZ-PP7 construct from Liu et al)

%% Variable inputs
if ~isempty(varargin)
    varargin = varargin{1};
end

%By default, don't use modified rate fitting
polyrate = false;
polyorder = 0;
meanrate = false;
prior_sigma = 0;
ratestring = '';
modelstring = '';
construct = 'P2P-MS2v5-LacZ-PP7v4'; %Base construct by default

for i=1:length(varargin)
    if strcmpi(varargin{i},'PolyRate')
        polyrate = true;
        polyorder = varargin{i+1};
        ratestring = ['_PolyOrder(n=',num2str(polyorder),')'];
    end
    if strcmpi(varargin{i},'MeanRate')
        meanrate = true;
        prior_sigma = varargin{i+1};
        ratestring = ['_MeanRate(s=',num2str(prior_sigma),')'];
    end
    if strcmpi(varargin{i},'ModelType')
        if strcmpi(varargin{i+1},'Dwell')
            modelstring = 'dwell';
            modeltypestring = 'dwell';
            modelloc = 'DwellModel\';
        elseif strcmpi(varargin{i+1},'Termination')
            modelstring = 'termrate';
            modeltypestring = 'termination';
            modelloc = 'TerminationModel\';
        end
    end
    if strcmpi(varargin{i},'Construct')
        construct = varargin{i+1};
    end
end

%Analysis type
if strcmp(analysistype,'initial')
    typestring = '_InitialRise';
elseif strcmp(analysistype,'whole')
    typestring = '_WholeCycle';
elseif strcmp(analysistype,'whole_using_initial')
    typestring = '_WholeUsingInitial';
end


%% Load data
%If using InitialRise data as given, choose which results to load.
%Otherwise, choose which dataset to infer results for.
initialresults_all = [];
if strcmp(analysistype,'whole_using_initial')
    initialresults = dir(['S:\Jonathan\Dropbox\2 Color Elongation\MCMC Results\',...
        modelloc,'*.mat']);
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
        initialresults_temp = load([initialresults(s(i)).folder,'\',initialresults_names{s(i)}]);
        initialresults_all(i).Prefix = initialresults_temp.Prefix;
        for j = 1:length(initialresults_temp.MCMCresults)
            initialresults_all(i).MCMCresults(j).mean_v = initialresults_temp.MCMCresults(j).mean_v;
            initialresults_all(i).MCMCresults(j).nucleus_index = initialresults_temp.MCMCresults(j).nucleus_index;
            initialresults_all(i).MCMCresults(j).ApprovedFits = initialresults_temp.MCMCresults(j).ApprovedFits;
            
        end
    end
    
    %Load the data
    fileloc = 'S:\Jonathan\Dropbox\2 Color Elongation\Filtered Particles';
    data_all = struct('ElapsedTime',{},'Prefix',{},'TwoColorFilteredParticles',{},...
        'nc12',{},'nc13',{},'nc14',{});
    for i = 1:length(initialresults_all)
        dat = load([fileloc,'\',initialresults_all(i).Prefix,'.mat']);
        temp.ElapsedTime = dat.ElapsedTime;
        temp.Prefix = dat.Prefix;
        temp.nc12 = dat.nc12;
        temp.nc13 = dat.nc13;
        temp.nc14 = dat.nc14;
        temp.TwoColorFilteredParticles = dat.TwoColorFilteredParticles;
        data_all(i) = temp;
    end
else
    files = dir('S:\Jonathan\Dropbox\2 Color Elongation\Filtered Particles\*.mat');
    names = {files.name};

    [s,v] = listdlg('PromptString','Select a dataset:','SelectionMode',...
        'multiple','ListString',names); %s is the indices of the chosen datasets
    
    %Load the data
    data_all = struct('ElapsedTime',{},'Prefix',{},'TwoColorFilteredParticles',{},...
        'nc12',{},'nc13',{},'nc14',{});
    for i = 1:length(s)
        dat = load([files(s(i)).folder,'\',names{s(i)}]);
        temp.ElapsedTime = dat.ElapsedTime;
        temp.Prefix = dat.Prefix;
        temp.nc12 = dat.nc12;
        temp.nc13 = dat.nc13;
        temp.nc14 = dat.nc14;
        temp.TwoColorFilteredParticles = dat.TwoColorFilteredParticles;
        data_all(i) = temp;
    end
end

w = waitbar(0,'Analyzing datasets...');

%% Load each dataset and set definitions
for k = 1:length(data_all)
waitbar(k/length(data_all),w,['Analyzing dataset ',num2str(k),' of ',num2str(length(data_all))]);

%Definitions
nc13flag = false;
nc13 = data_all(k).nc13;
nc14 = data_all(k).nc14;
t = data_all(k).ElapsedTime;
dt = mean(t((nc14+1):end)-t(nc14:(end-1)));
scalefac = 10; %To approximate # of RNAP's, divide fluorescence by this factor
datatype = '2-color';

%nc13 flag if nc13 mitosis wasn't recorded
if nc13 == 0
    nc13 = 1; %Set "fake" nc13 to frame 1
    nc13flag = true;
end

if strcmp(ncstring, 'nc14')
    nc = nc14;
elseif strcmp(ncstring, 'nc13')
    nc = nc13;
end

N = length(data_all(k).TwoColorFilteredParticles{1}); %Number of nuclei

%Setup data to be saved in structure
MCMCchain = struct('v_chain',{},'ton_chain',{},'A_chain',{},[modelstring,'_chain'],{},...
    'MS2_basal_chain',{},'PP7_basal_chain',{},'R_chain',{},'s2chain',{});
MCMCresults = struct('mean_v',{},'sigma_v',{},'mean_ton',{},'sigma_ton',{},...
    'mean_A',{},'sigma_A',{},['mean_',modelstring],{},['sigma_',modelstring],{},...
    'mean_MS2_basal',{},'sigma_MS2_basal',{},'mean_PP7_basal',{},'sigma_PP7_basal',{},'mean_R',{},...
    'sigma_R',{},'mean_sigma',{},'sigma_sigma',{},'nucleus_index',{},'ApprovedFits',{});
MCMCplot = struct('t_plot',{},'MS2_plot',{},'PP7_plot',{},'t_interp',{},...
    'MS2_interp',{},'PP7_interp',{},'simMS2',{},'simPP7',{},'MeanAP',{});
Prefix = data_all(k).Prefix;

%% MCMC analysis loop through single nuclei
parfor nuc = 1:N
PP7 = data_all(k).TwoColorFilteredParticles{1}(nuc).FilteredFluo./scalefac;
MS2 = data_all(k).TwoColorFilteredParticles{2}(nuc).FilteredFluo./scalefac;

%Check to make sure this particle exists in this nuclear cycle
if nc == nc13
    if sum(~isnan(MS2(nc13:nc14))) < 10
        disp(['skipping particle ',num2str(nuc), ' for nc13 analysis']);
        continue
    end
elseif nc == nc14
    if sum(~isnan(MS2(nc14:end))) < 10
        disp(['skipping particle ',num2str(nuc), ' for nc14 analysis']);
        continue
    end
end

%Define the analysis window: look at subset of nc13 or nc4
firstnc13time = 1.5; %Start 1.5 min after nc13
firstnc13timealt = 14.5; %If no nc13 captured, start 14.5 min before start of nc14
firstnc14time = 1.5; %Start 1.5 min after nc14

%Choose initial rise or whole cycle analysis
if strcmp(analysistype,'initial')
    lastnc13time = 8; %End 8 min before nc14
    lastnc14time = 10; %End 8 min after nc14
elseif strcmp(analysistype,'whole') || strcmp(analysistype,'whole_using_initial')
    lastnc13time = 5; %End 5 min before nc14
    lastnc14time = 18; %End 18 min after nc14
end

%Find the frames for interpolation cutoff times (separate for each signal). For the
%first timepoint, find the time of the first non-nan datapoint. For the last
%timepoint, find the time of the last non-nan datapoint.
if nc == nc13
    %First frame of interpolation
    if ~nc13flag
        firstnc13MS2 = find(and(t - t(nc13) > firstnc13time, ~isnan(MS2)),1);
        firstnc13PP7 = find(and(t - t(nc13) > firstnc13time, ~isnan(PP7)),1);
    else
        firstnc13MS2 = find(and(t(nc14) - t < firstnc13timealt, ~isnan(MS2)),1);
        firstnc13PP7 = find(and(t(nc14) - t < firstnc13timealt, ~isnan(PP7)),1);
    end
    %Last frame of interpolation
    lastnc13 = find(t(nc14) - t < lastnc13time,1);
    lastnc13MS2 = find(~isnan(MS2(1:lastnc13)),1,'last');
    lastnc13PP7 = find(~isnan(PP7(1:lastnc13)),1,'last');
elseif nc == nc14
    %First frame of interpolation
    firstnc14MS2 = find(and(t - t(nc14) > firstnc14time, ~isnan(MS2)),1);
    firstnc14PP7 = find(and(t - t(nc14) > firstnc14time, ~isnan(PP7)),1);
    %Last frame of interpolation
    lastnc14 = find(t - t(nc14) > lastnc14time,1);
    lastnc14MS2 = find(~isnan(MS2(1:lastnc14)),1,'last');
    lastnc14PP7 = find(~isnan(PP7(1:lastnc14)),1,'last');
end

%Now, interpolate the data.
if nc == nc13
    firstMS2 = firstnc13MS2;
    firstPP7 = firstnc13PP7;
    lastMS2 = lastnc13MS2;
    lastPP7 = lastnc13PP7;
    last = lastnc13;
elseif nc == nc14
    firstMS2 = firstnc14MS2;
    firstPP7 = firstnc14PP7;
    lastMS2 = lastnc14MS2;
    lastPP7 = lastnc14PP7;
    last = lastnc14;
end
t_interp = (t(nc):dt:t(last))-t(nc);

%Fill in nans with pchip interpolation
PP7_interp = fillmissing(PP7,'pchip');
MS2_interp = fillmissing(MS2,'pchip');

%Interpolate traces
PP7_interp = interp1(t(nc:last)-t(nc),...
    PP7_interp(nc:last),t_interp);
MS2_interp = interp1(t(nc:last)-t(nc),...
    MS2_interp(nc:last),t_interp);

%Define the subset of data to be analyzed. Analyze between the first
%and last datapoints in the nuclear cycle,
%but for the beginning cutoff leave the earlier times in, just with nan
%values. This is so we can still fit t_on, but we don't interpolate past
%the first datapoint.
if ~isempty(firstMS2)
    MS2_interp((t_interp+t(nc))<t(firstMS2)) = nan;
end
if ~isempty(firstPP7)
    PP7_interp((t_interp+t(nc))<t(firstPP7)) = nan;
end
if ~isempty(lastMS2)
    MS2_interp((t_interp+t(nc))>t(lastMS2)) = nan;
end
if ~isempty(lastPP7)
    PP7_interp((t_interp+t(nc))>t(lastPP7)) = nan;
end

%MCMC data structure
data = [];
data.xdata = t_interp; %Time data
if strcmp(datatype,'2-color')
    data.ydata = [MS2_interp,PP7_interp];
elseif strcmp(datatype,'1-color')
    data.ydata = MS2_interp;
end

%% Setup MCMC Fit

%Define the anonymous functions for the sum of squares function and model.
if polyrate
    ssfun = @(x,data) SumofSquaresConstant_MCMCstat_FreeScaling(construct,datatype,data,x,...
        {modeltypestring,x(4),x(1),'basal',x(5),x(6),'PolyRate'});
elseif meanrate
    ssfun = @(x,data) SumofSquaresConstant_MCMCstat_FreeScaling(construct,datatype,data,x,...
        {modeltypestring,x(4),x(1),'basal',x(5),x(6),'MeanRate'});
else
    ssfun = @(x,data) SumofSquaresConstant_MCMCstat_FreeScaling(construct,datatype,data,x,...
        {modeltypestring,x(4),x(1),'basal',x(5),x(6)});
end

% Initialize MCMC parameters

%Load previously inferred elongation rate if desired
if strcmp(analysistype,'whole_using_initial')
    nuctoload = find([initialresults_all(k).MCMCresults.nucleus_index] == nuc);
    v0 = initialresults_all(k).MCMCresults(nuctoload).mean_v;
    if isempty(v0)
        continue
    end
else
    v0 = 1+2*rand;
end
ton0 = 4*rand;
A0 = rand;
dwell0 = 4*rand;
termrate0 = 4*rand;
MS2_basal0 = 10;
PP7_basal0 = 5;

if polyrate
    R0 = 10*rand(1,polyorder); %Polynomial approximation coefficients
elseif meanrate
    R0 = [15, normrnd(0,3,1,size(MS2_interp,2)-1)]; %Mean rate + fluctuations
else
    R0 = 15+normrnd(0,3,1,size(MS2_interp,2)-1); %Non-approximated rate
end

if strcmpi(modeltypestring,'dwell')
    x0 = [v0,ton0,A0,dwell0,MS2_basal0,PP7_basal0,R0];
elseif strcmpi(modeltypestring,'termination')
    x0 = [v0,ton0,A0,termrate0,MS2_basal0,PP7_basal0,R0];
else
    warning('Please specify model type as Name/Value pair varargin');
end

sigma2_0 = 1; %Initial guess for error variance

%This is the variance in the parameter proposal distribution. Change these
%to change the proposal acceptance rate, or if the convergence doesn't look
%good.
if strcmp(analysistype,'whole_using_initial')
    v_step = 0.0000001;
else
    v_step = 0.05;
end
ton_step = dt;
A_step = 0.05;
dwell_step = 0.1;
termrate_step = 0.1;
MS2_basal_step = 1;
PP7_basal_step = 1;
R_step = 0.5*ones(size(R0));

if strcmpi(modeltypestring,'dwell')
    J0 = diag([v_step, ton_step, A_step, dwell_step, MS2_basal_step,...
    PP7_basal_step,R_step]); %Initial covariance matrix
elseif strcmpi(modeltypestring,'termination')
    J0 = diag([v_step, ton_step, A_step, termrate_step, MS2_basal_step,...
    PP7_basal_step,R_step]); %Initial covariance matrix
else
    error('Please specify model type as Name/Value pair varargin');
end

%% Setup MCMC parameters and options
%Limits on elongation rate (depending on if we're fixing it from previously
%inferred results or letting it be a free parameter)
if strcmp(analysistype,'whole_using_initial')
    v_lower = v0-0.00001;
    v_upper = v0+0.00001;
else
    v_lower = 0;
    v_upper = 10;
end
params = {
    {'v', x0(1), v_lower, v_upper}
    {'ton', x0(2), 0, 10}
    {'A', x0(3), 0, 1}
    {modelstring, x0(4), 0, 20}
    {'MS2_basal', x0(5), 0, 50}
    {'PP7_basal', x0(6), 0, 50}
    };

for i = 1:size(R0,2)
    if polyrate
        params{end+1} = {strcat('R',num2str(i)), x0(6+i), -100, 100};
    elseif meanrate
        if i == 1
            params{end+1} = {strcat('R',num2str(i)), x0(6+i), 0, 40}; %Mean rate
        else
            params{end+1} = {strcat('R',num2str(i)), x0(6+i), -30, 30, 0, prior_sigma}; %Fluctuations with priors around zero
        end
    else
        params{end+1} = {strcat('R',num2str(i)), x0(6+i), 0, 40};
    end
end

model = [];
model.ssfun = ssfun;
model.sigma2 = sigma2_0;
model.N = length(data.ydata);

%MCMC options
options = [];
options.nsimu = n_steps; %Number of steps
options.updatesigma = 1; %Update error variance
options.qcov = J0; %Initial covariance
options.burnintime = n_burn; %Burn in time
options.adaptint = 100;
options.method = 'dram';
options.verbosity = 0; %Decrease text output

%Run the MCMC
[results,chain,s2chain] = mcmcrun(model,data,params,options);

%Extract chain results into individual parameters
v_chain = chain(n_burn:end,1);
ton_chain = chain(n_burn:end,2);
A_chain = chain(n_burn:end,3);
modeltype_chain = chain(n_burn:end,4);
MS2_basal_chain = chain(n_burn:end,5);
PP7_basal_chain = chain(n_burn:end,6);
R_chain = chain(n_burn:end,7:end);

%Mean results
mean_v = mean(v_chain);
sigma_v = std(v_chain,1);
mean_ton = mean(ton_chain);
sigma_ton = std(ton_chain,1);
mean_A = mean(A_chain);
sigma_A = std(A_chain,1);
mean_modeltype = mean(modeltype_chain);
sigma_modeltype = std(modeltype_chain,1);
mean_MS2_basal = mean(MS2_basal_chain);
sigma_MS2_basal = std(MS2_basal_chain,1);
mean_PP7_basal = mean(PP7_basal_chain);
sigma_PP7_basal = std(PP7_basal_chain,1);
mean_R = mean(R_chain,1);
sigma_R = std(R_chain,1);
mean_sigma = sqrt(mean(s2chain));
sigma_sigma = std(sqrt(s2chain),1);

%Plotting variables
if polyrate
    [simMS2,simPP7] = GetFluorFromPolPos(construct,ConstantElongationSim(mean_v,mean_ton,...
        PolynomialRate(mean_R,data.xdata),data.xdata),...
    {modeltypestring,mean_modeltype,mean_v,'basal',mean_MS2_basal,mean_PP7_basal});
elseif meanrate
    [simMS2,simPP7] = GetFluorFromPolPos(construct,ConstantElongationSim(mean_v,mean_ton,...
        MeanRate(mean_R),data.xdata),...
    {modeltypestring,mean_modeltype,mean_v,'basal',mean_MS2_basal,mean_PP7_basal});
else
    [simMS2,simPP7] = GetFluorFromPolPos(construct,ConstantElongationSim(mean_v,mean_ton,mean_R,data.xdata),...
        {modeltypestring,mean_modeltype,mean_v,'basal',mean_MS2_basal,mean_PP7_basal});
end

simMS2 = mean_A * simMS2;
if nc == nc13
    t_plot = t(nc13:nc14) - t(nc13);
    MS2_plot = MS2(nc13:nc14);
    PP7_plot = PP7(nc13:nc14);
elseif nc == nc14
    t_plot = t(nc14:end) - t(nc14);
    MS2_plot = MS2(nc14:end);
    PP7_plot = PP7(nc14:end);
end


%{
%% Plot results

%Acceptance fraction
disp(strcat('Accepted proposal fraction was: ', num2str(1-results.rejected)));
close all;


%Chain visualization
figure();
mcmcplot(chain(:,1:6),[],results,'chainpanel');

%Histogram plots

%Elongation rate
figure();
mcmcplot(v_chain,[],[],'hist');
title('Elongation rate');

%Initiation time
figure();
mcmcplot(ton_chain,[],[],'hist');
title('Time of initiation onset');

%Calibration factor
figure();
mcmcplot(A_chain,[],[],'hist');
title('Calibration factor (ratio of MS2 / PP7 fluorescence)');

%Dwell time
figure();
mcmcplot(dwell_chain,[],[],'hist');
title('Dwell time (min)');

%Error variance
figure();
mcmcplot(sqrt(s2chain),[],[],'hist')
title('Error std posterior')

%Plot autocorrelations
figure();
plot(lags,auorr,'.-',lags([1 end]),[0 0],'k');
grid on
xlabel('lags')
ylabel('autocorrelation');
text(lags(end),0,sprintf('Effective Sample Size (ESS): %.0f_ ',ceil(mean(ESS))),'verticalalignment','bottom','horizontalalignment','right')
title('Markov Chain Auto Correlation')

%Plot loading rate inference
figure()
box on
hold on
fillyy(t_interp(1:end-1),mean_R+sigma_R,mean_R-sigma_R,[1 0.9 0.9]); %Error bars
%plot(t(1:end-1),R0,'b--'); %Initial guess
plot(t_interp(1:end-1),mean_R,'r-'); %Inferred fit
line([mean_ton,mean_ton],ylim,'Color','blue','LineWidth',1,'LineStyle','--');
legend('Inference error','Inferred fit');
title('Loading rate inference');

%Predictive envelope plot of model

figure();
hold on
h(1) = plot(t_plot,MS2_plot,'ro','DisplayName','MS2');
h(2) = plot(t_plot,PP7_plot,'go','DisplayName','PP7');
h(3) = plot(t_interp,MS2_interp,'r.-');
h(4) = plot(t_interp,PP7_interp,'g.-');
if strcmp(envelope,'slow')
    out = mcmcpred(results,chain,[],t_interp',modelfunMS2);
    mcmcpredplot(out);
    hold on
    out = mcmcpred(results,chain,[],t_interp',modelfunPP7);
    mcmcpredplot(out);
    legend('MS2','PP7',...
    '99%','95%','90%','50%','Inferred fit');
elseif strcmp(envelope,'fast')
    plot(t_interp,simMS2,'k-','LineWidth',1,'DisplayName','Inferred fit');
    plot(t_interp,simPP7,'k-','LineWidth',1);
    legend;
end
xlim([t_interp(1)*0.8, t_interp(end)*1.2]);
line([mean_ton,mean_ton],ylim,'Color','blue','LineWidth',1,'LineStyle','--'); %Time of initiation
line(xlim,mean_A*[mean_MS2_basal,mean_MS2_basal],'Color','red','LineWidth',1,'LineStyle','--'); %Basal MS2 fluorescence
line(xlim,[mean_PP7_basal,mean_PP7_basal],'Color','green','LineWidth',1,'LineStyle','--'); %Basal PP7 fluorescence
xlabel('Time'); ylabel('Fluorescence');
hold off
title('Prediction of the model')
uistack(h,'top');
%}


%{
%Corner plot
mcmc_chain_oneR(:,1,:) = [mcmc_chain(:,1), mcmc_chain(:,3),mcmc_chain(:,end-1)];
mcmc_chain_oneR = permute(mcmc_chain_oneR,[3 2 1]);
figure
ecornerplot(mcmc_chain_oneR,'ks',true,'color',[.6 .35 .3])
%}

%{
%% Text results
disp('MCMC Results');
disp('name      initial guess      mean inference     inference error');
disp(['Elongation rate (kb/min): ',num2str([v0, mean_v, sigma_v])]);
disp(['Time of initiation (min): ',num2str([ton0, mean_ton, sigma_ton])]);
disp(['Calibration factor (MS2/PP7): ',num2str([A0, mean_A, sigma_A])]);
disp(['Dwell time (min): ',num2str([dwell0, mean_dwell, sigma_dwell])]);
disp(['Basal MS2 fluorescence: ',num2str([MS2_basal0, mean_MS2_basal, sigma_MS2_basal])]);
disp(['Basal PP7 fluorescence: ',num2str([PP7_basal0, mean_PP7_basal, sigma_PP7_basal])]);
disp(['Measurement error: ',num2str([sigma2_0, mean_sigma, sigma_sigma])]);
%}

%% Save single nucleus data
MCMCchain(nuc).v_chain = v_chain;
MCMCchain(nuc).ton_chain = ton_chain;
MCMCchain(nuc).A_chain = A_chain;
MCMCchain(nuc).MS2_basal_chain = MS2_basal_chain;
MCMCchain(nuc).PP7_basal_chain = PP7_basal_chain;
MCMCchain(nuc).R_chain = R_chain;
MCMCchain(nuc).s2chain = s2chain;

MCMCresults(nuc).mean_v = mean_v;
MCMCresults(nuc).sigma_v = sigma_v;
MCMCresults(nuc).mean_ton = mean_ton;
MCMCresults(nuc).sigma_ton = sigma_ton;
MCMCresults(nuc).mean_A = mean_A;
MCMCresults(nuc).sigma_A = sigma_A;
MCMCresults(nuc).mean_MS2_basal = mean_MS2_basal;
MCMCresults(nuc).sigma_MS2_basal = sigma_MS2_basal;
MCMCresults(nuc).mean_PP7_basal = mean_PP7_basal;
MCMCresults(nuc).sigma_PP7_basal = sigma_PP7_basal;
MCMCresults(nuc).mean_R = mean_R;
MCMCresults(nuc).sigma_R = sigma_R;
MCMCresults(nuc).mean_sigma = mean_sigma;
MCMCresults(nuc).sigma_sigma = sigma_sigma;
MCMCresults(nuc).nucleus_index = nuc;

%If using previous results, carry over approval/rejection
if strcmp(analysistype,'whole_using_initial')
    MCMCresults(nuc).ApprovedFits = initialresults_all(k).MCMCresults(nuctoload).ApprovedFits;
else
    MCMCresults(nuc).ApprovedFits = 0; %No approval/rejection by default
end

%Model parameters
if strcmpi(modeltypestring,'dwell')
    MCMCchain(nuc).dwell_chain = modeltype_chain;
    MCMCresults(nuc).mean_dwell = mean_modeltype;
    MCMCresults(nuc).sigma_dwell = sigma_modeltype;
elseif strcmpi(modeltypestring,'termination')
    MCMCchain(nuc).termrate_chain = modeltype_chain;
    MCMCresults(nuc).mean_termrate = mean_modeltype;
    MCMCresults(nuc).sigma_termrate = sigma_modeltype;
end

MCMCplot(nuc).t_plot = t_plot;
MCMCplot(nuc).MS2_plot = MS2_plot;
MCMCplot(nuc).PP7_plot = PP7_plot;
MCMCplot(nuc).t_interp = t_interp;
MCMCplot(nuc).MS2_interp = MS2_interp;
MCMCplot(nuc).PP7_interp = PP7_interp;
MCMCplot(nuc).simMS2 = simMS2;
MCMCplot(nuc).simPP7 = simPP7;
MCMCplot(nuc).firstnc13time = firstnc13time;
MCMCplot(nuc).firstnc13timealt = firstnc13timealt;
MCMCplot(nuc).lastnc13time = lastnc13time;
MCMCplot(nuc).firstnc14time = firstnc14time;
MCMCplot(nuc).lastnc14time = lastnc14time;
MCMCplot(nuc).MeanAP = data_all(k).TwoColorFilteredParticles{1}(nuc).MeanAP;
end

%% Postprocessing: reject empty particle results, save dataset info
remove_indices = false(1,length(MCMCchain));
for nuc = 1:length(MCMCchain)
    if isempty(MCMCchain(nuc).v_chain)
        remove_indices(nuc) = true;
    end
end

MCMCchain(remove_indices) = [];
MCMCresults(remove_indices) = [];
MCMCplot(remove_indices) = [];

%% Save data into .mat structure
loc = ['S:\Jonathan\Dropbox\2 Color Elongation\MCMC Results\',modelloc];
filename = [date,'-',data_all(k).Prefix,'_',ncstring,ratestring,typestring];
save([loc,filename,'.mat'],'MCMCresults','MCMCplot','Prefix');

loc_rawchain = ['S:\Jonathan\Dropbox\2 Color Elongation\MCMC Results\',modelloc,'Raw Chains\'];
filename = [date,'-',data_all(k).Prefix,'_',ncstring,ratestring,typestring,'_RawChain'];
save([loc_rawchain,filename,'.mat'],'MCMCchain');
end

close(w);
disp(['MCMC analysis complete. Information stored in: ',loc]);
end

function R = PolynomialRate(a,t)
%Generates a polynomial rate function given coefficients a and time series
%t, of length size(t) - 1.
polyorder = length(a)-1;
R = zeros(size(t));
for i = 0:polyorder
    R = R + a(i+1) * t.^i;
end

%Remove last element of R to keep dimensions consistent in model
R = R(1:end-1);
end

function R = MeanRate(a)
%Generates a rate function given a mean rate and fluctuations.
%The mean rate is the first element of the input
%a, and the fluctuations are the other elements.
R = a(1) + a(2:end);
end
