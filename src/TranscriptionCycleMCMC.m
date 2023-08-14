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
%   'BayesianCoverage': set Bayesian coverage of the credible intervals of
%   the fitted parameter. Either set to value between 0 and 1, or use
%   '1std', '2std' or '3std' to set to a multiple of the coverage of the
%   standard deviation of a Gaussian
%   'loadMCMCsetup', true/false: option to import refined initiation values
%   and proposal parameters from previous MCMC results.
%   'loadPrevious', true/false: option to load previous inference results
%   to retain inferred elongation rate for hierarchical fit (see Liu et al,
%   Section S3.2
%   'construct', construct: option to specify a custom reporter gene that
%   must be defined in the subfunction library.m (default uses the
%   P2P-MS2-lacZ-PP7 construct from Liu et al)
%
%   Not yet implemented:
%   'MonteCarlo', MonteCarlo: constructs can be defined to
%   contain particular sites, where any polymerase can stall and
%   prematurely terminate. If MonteCarlo>0, then premature termination is
%   modeled as a Bernoulli experiment for each polymerase and the
%   likelihood function is approximated via Monte Carlo integration over
%   the simulated fluorescence dynamics. It is recommendable to use
%   MonteCarlo > 1000. If MonteCarlo==0, then the premature termination
%   dynamics is simulated as a deterministic reduction of polymerases at
%   stalling sites. (default = 0) %MG TODO: Consider stochastic simulation
%   of the whole transcription dynamics (not just premature termination).

%% 1) Variable inputs
%Default settings
fileDir = pwd;
saveLoc = pwd;
numParPools = 8;
n_burn = 30000;
n_adapt = 10000;
n_steps = 100000;
ratePriorWidth = 50;
PriorTrunc = 30;
t_start = 0;
t_end = Inf;
BayesianCoverage = 0.95; % Set Bayesian coverage of credible intervals
loadMCMCsetup = false;
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
        n_adapt = varargin{i+1};
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
    if strcmpi(varargin{i},'BayesianCoverage')
        if isa(varargin{i+1},'double')
            if and(varargin{i+1}>0,varargin{i+1}<1)
                BayesianCoverage = varargin{i+1};
            end
        elseif isa(varargin{i+1},'char')
            switch varargin{i+1}
                case '1std'
                    BayesianCoverage = 0.6627;
                case '2std'
                    BayesianCoverage = 0.9545;
                case '3std'
                    BayesianCoverage = 0.9973;
                otherwise
                    disp('Wrong input for Bayesian coverage. Set to default: 0.95');
            end
        else
            disp('Wrong input for Bayesian coverage. Set to default: 0.95');
        end
    end
    if strcmpi(varargin{i},'loadMCMCsetup')
        loadMCMCsetup = true;
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

%Error messages
if n_adapt>n_burn
    error('Too many adaptive steps. Use a value smaller than the number of burn-in steps.');
end

%% 2.1) Load construct details

%Query to the construct library regarding the construct details:
%'segments' and 'velocities' contain details about the elongation rates on
%different segments of the construct and the length of the construct;
%'stemloops' contain details about stem loop positions and contitutively bound fluorophores;
%'x_stall' contains the positions of stalling sites, where premature termination may happen upon stalling;
[ElongationSegments,stemloops,x_stall] = library(construct);

%Extract distinguishable names of velocity parameters
velocity_names = unique(ElongationSegments.velocities);

%% 2.2) Load data
data_all = LoadData(fileDir,loadPrevious,loadMCMCsetup,velocity_names,x_stall);

%% 3) Data analysis and export

%Set the likelihood function (observation model) as an anonymous function
switch MonteCarlo
    case 0
        %Set log-likelihood (sum-of-squares residual function)
        ssfun = @(x,data) getFluorescenceDynamicsSS(data,x,ElongationSegments,velocity_names,stemloops,x_stall);
    otherwise
        %Set Monte Carlo approximation of the likelihood
        error('Error: Monte Carlo option not yet implemented.');
        %ymodel = @(x,data) Likelihood(data,x,fixed_params) %MG: TODO
        %simulator = @(x,data) getFluorescenceDynamicsSS(data,x,ElongationSegments,velocity_names,stemloops,x_stall);
        %Or simulate MonteCarlo average: simulator = @(x,data)MonteCarloAverage(data,x,fixed_params) %MG: TODO
end

%% 3.1) Load each dataset and define export structures and save empty export structures

%Loop through different datasets
w = waitbar(0,'Analyzing datasets...'); %Start waitbar, showing the progress in terms of datasets

% Start parallel cluster and create temporary mat-files for each worker (to
% reduce the working memory for long mcmc chains)
defaultCluster = parcluster(parallel.defaultClusterProfile); %Get default cluster profile
maxWorkers = min(numParPools,defaultCluster.NumWorkers); % Use less or equal to numParPool workers (depending on availibility of workers)
spmd (maxWorkers)
    TempMatfileName = [tempname(saveLoc),'.mat']; % each worker gets a unique filename for its temporary mat-file
    TempMatfile = matfile(TempMatfileName, 'Writable', true);
end
% Create a parallel pool constant from the composite TempMatfile. This
% allows to use the same matfile with the same worker without the need
% to repeatedly transfer it to the worker
TempMatfileConstant = parallel.pool.Constant(TempMatfile);

try
for setNum = 1:length(data_all)

    waitbar(setNum/length(data_all),w,['Analyzing dataset ',num2str(setNum),' of ',num2str(length(data_all))]); %Update waitbar

    %Save k-th dataset into variable to avoid unnecessary communication
    %overhead to the parfor (otherwise the whole structure array 'data_all' is sent to each worker, not just 'data_all(k)')
    Dataset = data_all(setNum);

    N = length(Dataset.data); %Number of cells in dataset setNum

    %Setup empty structures to save data in (for each dataset!)
    [fields_params,fields_MCMCchain,fields_MCMCresults,MCMCchain,MCMCresults,MCMCplot] = generateExportStructures(N,velocity_names,x_stall);

    %Get metadata
    DatasetName = Dataset.data(1).name; %Assuming all the cells come from the same dataset
    Metadata.DatasetName = DatasetName;
    Metadata.FittedConstruct = construct;
    Metadata.ratePriorWidth = ratePriorWidth;
    Metadata.PriorTrunc = PriorTrunc;
    Metadata.t_start = t_start;
    Metadata.t_end = t_end;
    Metadata.BayesianCoverageCI = BayesianCoverage;
    Metadata.AdaptiveSteps = n_adapt;
    Metadata.BurnInSteps = n_burn;
    %Add number of MonteCarlo samples if MonteCarlo-version was used
    if and(MonteCarlo>0,not(isempty(x_stall)))
        Metadata.MonteCarlo = MonteCarlo;
    end

    % Save empty export structure for a single dataset into .mat structure
    % (Set paths for results (independent of operating system) and save
    % MCMC results and plots)
    save_filename = fullfile(saveLoc,[date,'-',DatasetName,'.mat']);
    save(save_filename,'MCMCresults','MCMCplot','Metadata','-v7.3'); %create a Version 7.3 MAT-file to enable/accelerate saving and loading parts of variables via matfile objects.

    %MCMC raw chains
    save_filename_chain = fullfile(saveLoc,[date,'-',DatasetName,'_RawChain','.mat']);
    save(save_filename_chain,'MCMCchain','-v7.3'); %create a Version 7.3 MAT-file to enable/accelerate saving and loading parts of variables via matfile objects.

    %% 3.2) MCMC analysis: Loop through single cells (paralellized)

    % Define a temporary .mat duplicate of the empty export structure for each
    % parallel worker. Each worker will store mcmc results parallely into its own
    % duplicate. Finally the parallely computed results will be stored in
    % the export structures predefined above. Although this can increase
    % the calculation time, it enables the parallel calculation of long
    % MCMC chains. Otherwise the required working memory could be beyond capacity.
    % spmd is an execution block where the enclosed code runs on all the
    % available parallel workers.
    spmd
        % Seed the variables in the matfile object
        TempMatfile.MCMCresults = cell(1,N);
        TempMatfile.MCMCplot = cell(1,N);
        TempMatfile.MCMCchain = cell(1,N);
        TempMatfile.gotResult = false(1, N);
    end

    %Clear empty export structure variables
    clear MCMCchain MCMCresults MCMCplot Metadata;

    %Do MCMC for each cell individually
    parfor cellNum = 1:N

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
        if loadMCMCsetup
            [params,~,~] = setupMCMC(fields_params,velocity_names,loadPrevious,t_interp,ratePriorWidth,PriorTrunc,x_stall);
            N_fields = length(fields_params)-1;
            N_params = length(params);
            for idxParam = 1:N_fields-1
                params{idxParam}{2} = Dataset.setup.(fields_MCMCresults{idxParam}) ;
            end
            setup_dR = Dataset.setup.mean_dR;
            for idxParam = 1:N_params-N_fields
                params{N_fields+idxParam}{2} = setup_dR(idxParam) ;
            end
            J0 = Dataset.setup.qcov;
            sigma2_initial = Dataset.setup.mean_sigma;
        else
            [params,J0,sigma2_initial] = setupMCMC(fields_params,velocity_names,loadPrevious,t_interp,ratePriorWidth,PriorTrunc,x_stall);
        end

        %Set 'model' (for mcmcrun)
        model = [];
        model.N = length(data.ydata);

        %Set the likelihood function (observation model) as an anonymous function for the sum of squares function and model.
        switch MonteCarlo
            case 0
                model.ssfun = ssfun;
                model.sigma2 = sigma2_initial;
            otherwise
                error('Error: Monte Carlo option not yet implemented.')
                %model.modelfun = ymodel; %MG: TODO
        end

        %Set 'options' (for mcmcrun):
        % The burn-in time is the number of steps we assign for convergence
        % of the MCMC to stationarity. Afterwards the adaptation interval
        % follows. The length of the exported chain is n_step-n_burn.
        options = [];
        options.nsimu = n_steps; %Number of steps
        options.updatesigma = 1; %s2 is not considered as fixed hyperparameter, but an actual parameter. The option sets the prior of s2 to be the conjugate prior (inverse gamma distribution parametrized as a scaled inverse chi^2 distribution with degree of freedom N0=1 and initial value of the scale parameter sigma2_initial)
        options.qcov = J0; %Initial covariance matrix of the proposal
        options.burnintime = n_burn; %Burn in time
        options.adaptint = n_adapt;
        options.method = 'dram';
        options.verbosity = 0; %Decrease text output

        %Run the MCMC
        [mcmcResult,chain,s2chain] = mcmcrun(model,data,params,options);

        %% 3.2.3) Extract and save results of MCMC for each single cell of the dataset
        if ~isempty(chain(n_burn:end,1)) % Check whether mcmcrun worked. If not, keep empty export structure for respective cell

            %Create empty export structures for single cells as temporary
            %variables (accessible only within a single iteration of the parfor loop)
            tempOne = sign(cellNum); %Get a temporary 1
            [~,~,~,localMCMCchain,localMCMCresults,localMCMCplot] = generateExportStructures(tempOne,velocity_names,x_stall);

            %Save the chains, the means of the resulting individual chains
            %and equal-tailed 95%-credible intervals
            N_idx=length(fields_params)+1;
            mean_params = zeros(1,N_idx-2); %Generate vector of mean parameters
            lower_quantile = (1 - BayesianCoverage)/2; upper_quantile = 1 - lower_quantile; % Set quantiles
            CredibleIntervals = quantile(chain(n_burn:end,:),[lower_quantile,upper_quantile]); %Get credible intervals

            % Loop through individual parameters
            for idx=1:(N_idx-2)
                localMCMCchain.(fields_MCMCchain{idx}) = chain(n_burn:end,idx);%Extract chain after burn in
                mean_params(idx) = mean(chain(n_burn:end,idx)); %Extract mean of chain after burn in
                localMCMCresults.(fields_MCMCresults{idx}) = mean_params(idx); %Save posterior mean estimate of the parameter
                localMCMCresults.(fields_MCMCresults{N_idx+idx}) = CredibleIntervals(:,idx); %Extract and save credible interval of the parameter
            end

            localMCMCchain.dR_chain = chain(n_burn:end,length(fields_params):end);
            mean_dR = mean(chain(n_burn:end,length(fields_params):end));
            localMCMCresults.mean_dR = mean_dR;
            localMCMCresults.CI_dR = CredibleIntervals(:,length(fields_params):end); %Save credible intervals of initiation fluctuations
            mean_params = [mean_params,mean_dR]; %Add initiation fluctuations to mean parameter vector (-> vector required for simulation of best fit)

            localMCMCchain.s2chain = s2chain;
            localMCMCresults.mean_sigma = sqrt(mean(s2chain(n_burn:end))); %square root of posterior mean of error variance s2
            localMCMCresults.CI_sigma = sqrt(quantile(s2chain(n_burn:end),[lower_quantile,upper_quantile]))'; %square root of credible interval of error variance s2

            localMCMCresults.cell_index = cellNum; %cell number
            localMCMCresults.mcmcrun = mcmcResult; %data about mcmc


            %If using previous results, carry over approval/rejection
            if loadPrevious
                localMCMCresults.ApprovedFits = initialresults_all(setNum).MCMCresults(cellToload).ApprovedFits;
            else
                localMCMCresults.ApprovedFits = 0; %No approval/rejection by default
            end

            %Simulated fluorescences of best fit (with posterior mean parameters)
            [~,simMS2,simPP7] = getFluorescenceDynamicsSS(data,mean_params,ElongationSegments,velocity_names,stemloops,x_stall);

            % Save data and best fit (from t_start to t_end) for plotting
            localMCMCplot.t_plot = t;
            localMCMCplot.MS2_plot = MS2;
            localMCMCplot.PP7_plot = PP7;
            localMCMCplot.simMS2 = simMS2;
            localMCMCplot.simPP7 = simPP7;

            % Save all chains, results and plots into temporary .mat files.
            % Afterwards the MCMC results of the particular cell are
            % discarded from the working memory. This is particularly
            % useful for a large number of long chains.
            matfileObj = TempMatfileConstant.Value; % Access the actual matfile object 'TempMatfile' associated with the worker, where the loop iteration is currently running on.
            matfileObj.MCMCchain(1,cellNum) = {localMCMCchain};
            matfileObj.MCMCresults(1,cellNum) = {localMCMCresults};
            matfileObj.MCMCplot(1,cellNum) = {localMCMCplot};
            matfileObj.gotResult(1,cellNum) = true; % Report that results for this particular cellNum have been stored into this particular temporary .mat file
        end
    end %End of (paralellized) loop through cells

    % Get matfile objects for the final export MAT-files
    m = matfile(save_filename, 'Writable', true);
    mCHAIN =  matfile(save_filename_chain, 'Writable', true);

    % Sequentially transfer results from temporary MAT-files to final
    % MAT-files (cell by cell)
    for worker_idx = 1: numel(TempMatfileName)
        workerTempMatfileName = TempMatfileName{worker_idx};
        workerTempMatfile = matfile(workerTempMatfileName);
        for cellNum = 1:N
            if workerTempMatfile.gotResult(1,cellNum)
                localCHAINS = workerTempMatfile.MCMCchain(1,cellNum);
                mCHAIN.MCMCchain(1,cellNum) = localCHAINS{1};
                localRESULTS = workerTempMatfile.MCMCresults(1,cellNum);
                m.MCMCresults(1,cellNum) = localRESULTS{1};
                localPLOT = workerTempMatfile.MCMCplot(1,cellNum);
                m.MCMCplot(1,cellNum) = localPLOT{1};
            end
        end
    end
end %End of loop through datasets
end %End of try

% Delete the temporary MAT-files after using them
for worker_idx = 1: numel(TempMatfileName)
    delete(TempMatfileName{worker_idx});
end
close(w); % Close waitbar

disp(['MCMC analysis complete. Information stored in: ',saveLoc]);
end
