function [params,J0,sigma2_initial] = setupMCMC(fields_params,velocity_names,loadPrevious,t,ratePriorWidth,PriorTrunc,x_stall)
% setupMCMC fixes the parameters for inference, its priors, initial values of the MCMC chain and the covariance matrix of the proposal distribution of the MCMC.
% params: input argument of mcmcrun, specifying the parameter names, its initial values and the marginal prior distributions
% J0: covariance matrix of the proposal distribution
% sigma2_initial: initial value of the observational parameter s2 in the MCMC chain
%
% ratePriorWidth,PriorTrunc: parameters of the truncated Gaussian prior for the parameters dR

%% Initial values of MCMC chains.
MCMC_initial = cell2struct(cell(1,length(fields_params)),fields_params,2);
% JL 10/28/2020: Need to update this to allow user specification of MCMC
% parameters, bounds, and priors

%Set initial values of the MCMC chains for the elongation rates. Load
%previously inferred elongation rate if desired.
if loadPrevious
    %MG 04/18/2023: Need to adapt this to allow for loading of previous
    %data (hereby currently disabled)
    %{
    cellToload = find([initialresults_all(k).MCMCresults.cell_index] == cellNum);
    v0 = initialresults_all(k).MCMCresults(cellToload).mean_v;

    if isempty(MCMC_initial.v)
        continue
    end
    %}
else
    for idx=1:length(velocity_names)
        MCMC_initial.(velocity_names{idx}) = 1+2*rand; %Initial guess for elongation rate (kb/min)
    end
end

%Initials for other transcription cycle parameters
MCMC_initial.ton = 4*rand;
MCMC_initial.A = rand;
MCMC_initial.tau = 4*rand;
MCMC_initial.MS2_basal = 10;
MCMC_initial.PP7_basal = 5;
MCMC_initial.R = 15;
MCMC_initial.dR = normrnd(0,3,1,length(t)); %The number of initiation fluctuation parameters equals the length of t_interp (sic!)

%Initials for premature termination parameters
if not(isempty(x_stall))
    MCMC_initial.ProbPremTerm = 0.001; %Choose initial value, such that low probabilities are favored
    %MCMC_initial.ProbPremTerm = 0.5; %Choose initial value to represent
    %maximum ignorance about probability of premature termination
    MCMC_initial.tauPremTerm = 4*rand;
end
sigma2_initial = 0.2; %Initial value of the MCMCchain for the error variance parameter s2 of the nonlinear Gaussian observation model

%% Set step size (covariance matrix) of the MCMC proposal distribution.
%Change these to change the proposal acceptance rate, or if the convergence
%doesn't look good.
MCMC_step = cell2struct(cell(1,length(fields_params)),fields_params,2);
if loadPrevious
    for idx=1:length(velocity_names)
        MCMC_step.(velocity_names{idx}) = 0.0000001;
    end
else
    for idx=1:length(velocity_names)
        MCMC_step.(velocity_names{idx}) = 0.05;
    end
end
MCMC_step.ton = t(end)-t(end-1);
MCMC_step.A = 0.05;
MCMC_step.tau = 0.1;
MCMC_step.MS2_basal = 1;
MCMC_step.PP7_basal = 1;
MCMC_step.R = 0.5;
MCMC_step.dR = 0.5*ones(size(MCMC_initial.dR));

%Step size for premature termination parameters
if not(isempty(x_stall))
    MCMC_step.ProbPremTerm = 0.05;
    MCMC_step.tauPremTerm = 0.1;
end

%Feed initial step sizes into a (covariance) matrix (requirement of mcmcrun)
proposalstep = zeros(1,length(fields_params)-1+length(MCMC_step.dR));
for idx=1:(length(fields_params)-1)
    proposalstep(idx) = MCMC_step.(fields_params{idx});
end
proposalstep(length(fields_params):end) = MCMC_step.dR;
J0 = diag(proposalstep); %Initial covariance matrix

%% Set 'params': specifications parameters of interest (parameter names, its initial values and the marginal prior distributions)

%Bounds of elongation rate uniform prior (depending on if we're fixing it from previously
%inferred results or letting it be a free parameter)

%Set lower and upper bounds for velocities
v_lower=zeros(1,length(velocity_names));
v_upper=ones(1,length(velocity_names));
if loadPrevious
    for idx=1:length(velocity_names)
        v_lower(idx) = MCMC_initial.(velocity_names{idx}) - 0.00001;
        v_upper(idx) = MCMC_initial.(velocity_names{idx}) + 0.00001;
    end
else
    %v_lower = v_lower * 0;
    v_upper = v_upper * 10; %elongation rate priors are homogeneous
end

%Set mcmcrun input argument "params" to specify parameter names, initial
%values and priors.
params = cell(length(fields_params)-1+length(MCMC_step.dR),1); % initialize cells array as a column

%Set elongation rate specifications
N_idx = length(velocity_names);
for idx=1:N_idx
    params{idx} = {fields_params{idx}, MCMC_initial.(velocity_names{idx}),v_lower(idx),v_upper(idx)}; %Add elongation rate parameters to params
end

%Set specifications for other transcription cycle parameters
params{N_idx+1} = {'tau', MCMC_initial.tau, 0, 20};
params{N_idx+2} = {'ton', MCMC_initial.ton, 0, 10};
params{N_idx+3} = {'MS2_basal', MCMC_initial.MS2_basal, 0, 50};
params{N_idx+4} = {'PP7_basal', MCMC_initial.PP7_basal, 0, 50};
params{N_idx+5} = {'A', MCMC_initial.A, 0, 1};
params{N_idx+6} = {'R', MCMC_initial.R, 0, 40};

%Add premature termination parameters to params
if not(isempty(x_stall))
    params{N_idx+7} = {'ProbPremTerm', MCMC_initial.ProbPremTerm, 0, 1};
    params{N_idx+8} = {'tauPremTerm', MCMC_initial.tauPremTerm, 0, 20};
end

%Set initial values and priors for initiation rate fluctuation parameters
N_idx = length(fields_params)-1;
for k = 1:size(MCMC_initial.dR,2)
    params{N_idx+k} = {strcat('dR',num2str(k)), MCMC_initial.dR(k), -PriorTrunc, PriorTrunc, 0, ratePriorWidth}; %Fluctuations with priors around zero
end