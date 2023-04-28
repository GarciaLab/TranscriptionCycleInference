function TranscriptionCycleMCMC(varargin)

%Fits transcription cycle model from Liu et al (2020) using Markov Chain
%Monte Carlo, implemented in the MCMCstat package. The data must be
%generated using a dual color MS2/PP7 reported as described in Liu et al.

%The code can load multiple datasets, and possesses user-specified options
%as variable arguments. The details of these inputs are described in the
%readme and in the manuscript, and are initialized with default values if
%not specified by the user.


%Variable arguments:
%   'fileDir', fileDir: directory of dataset files (default = root
%   directory)
%   'saveLoc', saveLoc: directory of saved MCMC results (default = root
%   directory)
%   'numParPools', numParPools: number of parallel workers to use (default
%   = 8)
%   'n_burn', n_burn: number of burn-in steps in MCMC algorithm (default = 10000)
%   'n_steps', n_steps: number of steps in MCMC algorithm including burn-in (default =
%   20000)
%   'ratePriorWidth', ratePriorWidth: standard deviation of (truncated) Gaussian
%   prior for rate fluctuation term dR(t) (default = 50, see Liu et al, Section S3.1)
%   'PriorTrunc', PriorTrunc: symmetric truncation of the prior for rate
%   fluctuation term dR(t) to the interval (-PriorTrunc, PriorTrunc). (default = 30)
%   't_start', t_start: time value to start fit at (default = 0)
%   't_end', t_end: time value to end fit at (default = Inf)
%   'loadPrevious', true/false: option to load previous inference results
%   to retain inferred elongation rate for hierarchical fit (see Liu et al,
%   Section S3.2
%   'construct', construct: option to specify a custom reporter gene that
%   must be defined in the subfunction GetFluorFromPolPos (default uses the
%   P2P-MS2-lacZ-PP7 construct from Liu et al)
%
%   'MonteCarlo', MonteCarlo: constructs can be defined to
%   contain particular sites, where any polymerase can drop off the
%   construct leading to premature termination. If MonteCarlo>0,
%   then drop-off is modeled as a Bernoulli experiment for each polymerase
%   and the likelihood function is approximated via Monte Carlo integration
%   over the simulated fluorescence dynamics. It is recommendable to use
%   MonteCarlo > 1000. If MonteCarlo==0, then the drop-off
%   dynamics is simulated as a deterministic reduction of polymerases at
%   drop-off sites. (default = 0) %MG TODO: Consider stochastic simulation
%   of the whole transcription dynamics (not just drop-off).

%% 1) Variable inputs
%Default settings
fileDir = pwd;
saveLoc = pwd;
numParPools = 8;
n_burn = 30000;
n_adapt = 200;
n_steps = 100000;
ratePriorWidth = 50;
PriorTrunc = 30;
t_start = 0;
t_end = Inf;
loadPrevious = false;
construct = 'P2P-MS2v5-LacZ-PP7v4';
MonteCarlo = 0;

for i=1:length(varargin)
    if strcmpi(varargin{i},'fileDir')
        fileDir = varargin{i+1};
    end
    if strcmpi(varargin{i},'saveLoc')
        saveLoc = varargin{i+1};
    end
    if strcmpi(varargin{i},'numParPools')
        numParPools = varargin{i+1};
    end
    if strcmpi(varargin{i},'n_burn')
        n_burn = varargin{i+1};
    end
    if strcmpi(varargin{i},'n_adapt')
        n_burn = varargin{i+1};
    end
    if strcmpi(varargin{i},'n_steps')
        n_steps = varargin{i+1};
    end
    if strcmpi(varargin{i},'ratePriorWidth')
        ratePriorWidth = varargin{i+1};
    end
    if strcmpi(varargin{i},'PriorTrunc')
        construct = varargin{i+1};
    end
    if strcmpi(varargin{i},'t_start')
        t_start = varargin{i+1};
    end
    if strcmpi(varargin{i},'t_end')
        t_end = varargin{i+1};
    end
    if strcmpi(varargin{i},'loadPrevious')
        loadPrevious = true;
    end
    if strcmpi(varargin{i},'construct')
        construct = varargin{i+1};
    end
    if strcmpi(varargin{i},'MonteCarlo')
        construct = varargin{i+1};
    end
end

%% 2.1) Load construct details

%Query to the construct library regarding the construct details:
%'segments' and 'velocities' contain details about the elongation rates on
%different segments of the construct and the length of the construct;
%'stemloops' contain details about stem loop positions and contitutively bound fluorophores;
%'x_drop' contains the positions of potential drop-off sites;
[ElongationSegments,stemloops,x_drop] = library(construct);

%Extract distinguishable names of velocity parameters
velocity_names = unique(ElongationSegments.velocities);

%% 2.2) Load data
data_all = LoadData(fileDir,loadPrevious);

%% 3) Data analysis and export

%Set the likelihood function (observation model) as an anonymous function
switch MonteCarlo
    case 0
        %Set log-likelihood (sum-of-squares residual function)
        ssfun = @(x,data) getFluorescenceDynamicsSS(data,x,ElongationSegments,velocity_names,stemloops,x_drop,true);
        simulator = @(data,x) getFluorescenceDynamicsSS(data,x,ElongationSegments,velocity_names,stemloops,x_drop,false);
    otherwise
        %Set Monte Carlo approximation of the likelihood
        error('Error: Drop off dynamics not yet implemented.');
        %ymodel = @(x,data) Likelihood(data,x,fixed_params) %MG: TODO
        %simulator = @(x,data) getFluorescenceDynamicsSS(data,x,ElongationSegments,velocity_names,stemloops,x_drop,false);
        %Or simulate MonteCarlo average: simulator = @(x,data)
        %MonteCarloAverage(data,x,fixed_params) %MG: TODO
end

%% 3.1) Load each dataset and define export structures

%Loop through different datasets
w = waitbar(0,'Analyzing datasets...'); %Start waitbar, showing the progress in terms of datasets
for k = 1:length(data_all)

    waitbar(k/length(data_all),w,['Analyzing dataset ',num2str(k),' of ',num2str(length(data_all))]); %Update waitbar

    %Save k-th dataset into variable to avoid unnecessary communication
    %overhead to the parfor (otherwise the whole structure array 'data_all' is sent to each worker, not just 'data_all(k)')
    Dataset = data_all(k);

    N = length(Dataset.data); %Number of cells in dataset k

    %Setup empty structures to save data in (for each dataset!)
    [fields_params,fields_MCMCchain,fields_MCMCresults,MCMCchain,MCMCresults,MCMCplot] = generateExportStructures(N,velocity_names,x_drop);

    %Get name of the dataset
    DatasetName = Dataset.data(1).name; %Assuming all the cells come from the same dataset

    %% 3.2) MCMC analysis: Loop through single cells (paralellized)

    %Do MCMC for each cell individually
    parfor (cellNum = 1:N, numParPools)

        %% 3.2.1) Load single cell data
        %Define time vector
        t = Dataset.data(cellNum).time; %Time vector for this cell

        %Initialize MS2 and PP7 signals
        MS2 = Dataset.data(cellNum).MS2;
        PP7 = Dataset.data(cellNum).PP7;

        %Truncate times to user-specified range
        indStart = find(t >= t_start,1,'first'); %Index of first used timepoint
        indEnd = find(t < t_end,1,'last'); %Index of last used timepoint

        t = t(indStart:indEnd); %Truncate t
        dt = mean(t(2:end)-t(1:(end-1))); %Average time resolution of the dataset for model usage
        t_interp = t(1):dt:t(end); %Interpolated time vector with even time resolution (for simulation)


        MS2 = MS2(indStart:indEnd);
        PP7 = PP7(indStart:indEnd);

        %MCMC 'data' structure initialization (input argument of mcmcrun)
        data = [];
        data.xdata = {t,t_interp}; %Time data and interpolated times
        data.ydata = [MS2,PP7]; %Fluorescence data

        %% 3.2.2) Setup MCMC Fit for single cell
        % The function mcmcrun requires four input arguments:
        % model: likelihood model and associated (hyper)parameters (number of observation points, initial value of s2)
        % data: already set in section 3.2.1)
        % params: parameter specifications (names of parameters, initial values of the MCMC chains and specification of the marginal priors)
        % options: settings for the MCMC (length of MCMC chain, burn in steps, adaptive method, covariance matrix of the proposal distribution, options for the observational parameter s2 set to sigma2_initial.)


        % Set 'params', the covariance matrix of the proposal distribution
        % 'J0' and the initial value 'sigma2_initial' of the observational parameter s2
        [params,J0,sigma2_initial] = setupMCMC(fields_params,velocity_names,loadPrevious,t_interp,ratePriorWidth,PriorTrunc,x_drop);

        %Set 'model' (for mcmcrun)
        model = [];
        model.N = length(data.ydata);

        %Set the likelihood function (observation model) as an anonymous function for the sum of squares function and model.
        switch MonteCarlo
            case 0
                model.ssfun = ssfun;
                model.sigma2 = sigma2_initial;
            otherwise
                error('Error: Drop off dynamics not yet implemented.')
                %model.modelfun = ymodel; %MG: TODO
        end

        %Set 'options' (for mcmcrun):
        % The burn-in time is the number of steps we assign for convergence
        % of the MCMC to stationarity. Afterwards the adaptation interval
        % follows. The length of the exported chain is n_step-n_burn.
        options = [];
        options.nsimu = n_steps; %Number of steps
        options.updatesigma = 1; %s2 is not considered as fixed hyperparameter, but an actual parameter. The option sets the prior of s2 to be the conjugate prior (inverse gamma distribution parametrized as a scaled inverse chi^2 distribution with degree of freedom N0=1 and initial value of the cale parameter sigma2_initial)
        options.qcov = J0; %Initial covariance matrix of the proposal
        options.burnintime = n_burn; %Burn in time
        options.adaptint = n_adapt;
        options.method = 'dram';
        options.verbosity = 0; %Decrease text output

        %Run the MCMC
        [~,chain,s2chain] = mcmcrun(model,data,params,options);

        %% 3.2.3) Extract and save results of MCMC for each single cell of the dataset
        % Save MCMC chains for individual parameters into the structure MCMCchain
        for idx=1:(length(fields_MCMCchain)-2)
            MCMCchain(cellNum).(fields_MCMCchain{idx}) = chain(n_burn:end,idx);
        end
        MCMCchain(cellNum).dR_chain = chain(n_burn:end,length(fields_params):end);
        MCMCchain(cellNum).s2chain = s2chain;

        %Save cell number and fitted construct
        MCMCresults(cellNum).cell_index = cellNum;
        MCMCresults(cellNum).FittedConstruct = construct;

        %Save means of the resulting individual chains and upper bound to
        %the error of the posterior mean estimator (Geyer's monotone sequence estimator)
        %MG TODO: so far the multivariate(!) monotone sequence estimator is not implemented and instead the standard deviation with normalization N is returned)
        N_idx=(length(fields_MCMCresults)-3)/2;
        mean_params = zeros(1,N_idx-2); %Generate vector of mean parameters
        for idx=1:(N_idx-2)
            mean_params(idx) = mean(MCMCchain(cellNum).(fields_MCMCchain{idx}));
            MCMCresults(cellNum).(fields_MCMCresults{idx}) = mean_params(idx); %Extract and save posterior mean of the parameter
        end
        MCMCresults(cellNum).mean_dR = mean(MCMCchain(cellNum).dR_chain);
        mean_params = [mean_params,MCMCresults(cellNum).mean_dR]; %Add initiation fluctuations to mean parameter vector
        MCMCresults(cellNum).mean_sigma = sqrt(mean(MCMCchain(cellNum).s2chain)); %square root of posterior mean of error variance s2

        for idx=1:(N_idx-1)
            MCMCresults(cellNum).(fields_MCMCresults{N_idx+idx}) = std(MCMCchain(cellNum).(fields_MCMCchain{idx}),1);
            %both std(...) and std(...,1) are taken along columns, but the
            %first normalizes the variance by (N-1) and the second by N.
        end
        MCMCresults(cellNum).sigma_sigma = sqrt(std(MCMCchain(cellNum).s2chain,1)); %square root of posterior standard deviation of error variance s2


        %If using previous results, carry over approval/rejection
        if loadPrevious
            MCMCresults(cellNum).ApprovedFits = initialresults_all(k).MCMCresults(cellToload).ApprovedFits;
        else
            MCMCresults(cellNum).ApprovedFits = 0; %No approval/rejection by default
        end

        %Simulated fluorescences of best fit (with posterior mean parameters)
        [simMS2,simPP7] = simulator(data,mean_params);

        % Save data and best fit (from t_start to t_end) for plotting
        MCMCplot(cellNum).t_plot = t;
        MCMCplot(cellNum).MS2_plot = MS2;
        MCMCplot(cellNum).PP7_plot = PP7;
        MCMCplot(cellNum).simMS2 = simMS2;
        MCMCplot(cellNum).simPP7 = simPP7;
    end %End of (paralellized) loop through cells

    %Add number of MonteCarlo samples if MonteCarlo-version was used
    if and(MonteCarlo>0,not(isempty(x_drop)))
        for idx=1:N
            MCMCresults(idx).MonteCarlo = MonteCarlo;
        end
    end

    %% 4.1) Postprocessing: reject empty particle results, save dataset info
    remove_indices = false(1,length(MCMCchain));
    for cellNum = 1:length(MCMCchain)
        if isempty(MCMCchain(cellNum).R_chain)
            remove_indices(cellNum) = true;
        end
    end

    MCMCchain(remove_indices) = [];
    MCMCresults(remove_indices) = [];
    MCMCplot(remove_indices) = [];

    %% 4.2) Save results for a single dataset into .mat structure
    %Set paths for results (independent of operating system) and save
    %MCMC results and plots
    filename = [date,'-',DatasetName,'.mat'];
    save(fullfile(saveLoc,filename),'MCMCresults','MCMCplot','DatasetName');

    %MCMC raw chains
    filename = [date,'-',DatasetName,'_RawChain','.mat'];
    save(fullfile(saveLoc,filename),'MCMCchain');
end %End of loop through datasets

close(w);

disp(['MCMC analysis complete. Information stored in: ',saveLoc]);

end