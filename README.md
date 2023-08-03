# TranscriptionCycleInference
Repository for transcription cycle MCMC inference code, developed on the basis of the code used in Liu et. al. (2020) (https://www.biorxiv.org/content/10.1101/2020.08.29.273474v3).

# Introduction
This repository contains the platform-independent MATLAB code for running the Markov Chain Monte Carlo (MCMC) fitting on dual-color fluorescent reporter data, as described in Liu et. al. (2020). This code utilizes the MCMCStat package by Haario. et al., available here: https://github.com/mjlaine/mcmcstat. It returns the posterior mean estimates (and equally-tailed 95% credible intervals) of the transcription cycle parameters specified for the reporter construct.

# Organization
The necessary scripts to run can be found in the \src folder. The primary script is `TranscriptionCycleMCMC.m`, which contains the master code that handles the MCMC fits. The file `library.m` contains (i) the information on the reporter constructs needed to calculate fluorescent signals from the transcription cycle model, (ii) the assignment of different elongation rate parameters to different segments of the construct, and (iii) the potential specification of stalling sites, where transcription may terminate prematurely. The user can add new reporter constructs to library.m to customize the inference. `ApproveMCMCResults.m` is a manual curation script that can be used to visualize fit results and approve/reject them. It includes the options of (i) visualizing histograms and MCMC chains of certain parameters and (ii) the fluctuations of the initiation rate. The \src\dependencies folder contains sub-functions that are necessary, but should not be altered by the user. The function `GetEvidence.m` can be used to obtain the evidence of a dataset with respect to a specific model (likelihood & prior), allowing model comparison via Bayes factor.

# MCMC Fitting Guide
For more details on the MCMC approach, refer to Section S3 of Liu et. al. (2020). To run an MCMC fit, TranscriptionCycleMCMC.m must be executed on a dataset, which must be saved as a MATLAB .mat file. The required data format for each dataset is as follows: there must be a structure array named <code>data</code>. Each index of <code>data</code> corresponds to a single cell of fluorescent data. Each index (i.e. single cell) of <code>data</code> must have four fields:

- <code>time</code>: a 1 x N vector of time points of experimental acquisition
- <code>MS2</code>: a 1 x N vector of MS2 fluorescent signals
- <code>PP7</code>: a 1 x N vector of PP7 fluorescent signals
- <code>name</code>: a string used to label this dataset
  
TranscriptionCycleMCMC.m loads one or more datasets with this format and runs the MCMC inference procedure on them. To customize the fitting, the following variable arguments can be passed:

- 'fileDir' (string): Directory to search for dataset files (default: current working directory).
- 'saveLoc' (string): Directory to save MCMC results (default: current working directory).
- 'numParPools' (int): Number of parallel pool workers (default: 8).
- 'n_burn' (int): Burn-in steps for MCMC (default: 30000).
- 'n_adapt' (int): Adaptive steps of adaptive MCMC (default: 10000)
- 'n_steps' (int): Total MCMC steps, including burn-in (default: 100000).
- 'ratePriorWidth' (float): Standard deviation of the (truncated) Gaussian prior for rate fluctuation terms dR(t) (default: 50, see Section S3.1).
- 'PriorTrunc' (int): Symmetric truncation of the prior for rate fluctuation terms dR(t) to the interval (-PriorTrunc, PriorTrunc) (default: 30).
- 't_start' (float): Time value to start fit at (default: 0)
- 't_end' (float): Time value to end fit at (default: Inf)
- 'BayesianCoverage' (float): Coverage of equally-tailed Bayesian credible intervals of the fitted parameters.
- 'loadPrevious' (boolean): Option to load previous inference results to retain inferred elongation rate for hierarchical fit (default: false, see Section S3.2)
- 'construct' (string): Specify a custom reporter gene defined in the subfunction library.m (default: 'P2P-MS2-lacZ-PP7' construct from Liu et. al. 2020)
<!-- - 'MonteCarlo' (int): option of approximate Bayesian computing via Monte Carlo approximation of a probabilistic elongation model -->

The results will be saved in two .mat files, one with the label 'RawChains' that contains the raw MCMC chains, and another that contains the summary statistics (posterior mean, equally-tailed credible intervals, etc.) and relevant plotting variables, contained in data structures named MCMCresults and MCMCplot, respectively.

# Defining Custom Reporter Constructs
To use a different reporter gene than the P2P-MS2-lacZ-PP7 construct in Liu et. al. 2020, simply define a construct in `library.m`. There is a template provided within `library.m`, where you can edit the following specifications:
- Positions of the MS2 and PP7 stem loop sequences as well as the number of stem loops per sequence.
- Number of fluorophores bound consistently (to the gene itself) throughout the observation period.
- Length of the construct and assignment of distinct elongation rate parameters to subsegments of the construct.
- Stalling sites, where a fraction of polymerases terminates prematurely. Setting stalling sites alters the model used by `TranscriptionCycleInference.m`. If you prefer to analyze stalling without incorporating premature termination into the likelihood model, consider setting distinct elongation rate parameters near the construct's stalling sites.

Ensure that your custom reporter gene has a unique label. Then pass this label into `TranscriptionCycleMCMC.m` as an argument using the 'construct' variable.

# Curating results: Approve/reject fits upon visual inspection
The function `ApproveMCMCResults.m` can be used to visualize MCMC fit results and optionally approve or reject fits manually. Running the script will open a dialog box to load the MCMC results data structure created from running `TranscriptionCycleMCMC.m`. The following variable arguments can be passed:

- 'fileDir' (string): directory to load MCMC results from (default: current working directory).
- 'RawChains': view raw MCMC samples in addition to mean fit results (default: no raw chains)
- 'RawChainsCheckbox': view raw MCMC chains and select the parameters via a dialog box (default: no raw chains)
- 'InitiationFluctuations': option to simultaneously view the initiation rate fluctuations (default: no initiation fluctioations)
- 'SmoothRate' (float): visualize time-dependent initiation rate using a user-defined smoothing span (e.g. span = 0.1 will smooth the initiation rate fluctuations over 10% of overall time window) (default = 0.1)
- 'LoadPrevious': import the results of a previous data curation (ApprovedFits) and overwrite the respective structures in the currently loaded dataset.

With this script, you can flip between different single-cell fits and visualize the data with mean fit results, as well as the raw MCMC chains and their histograms. You can then approve or reject a particular fit result manually, and save these choices in a field of the MCMCresults data structure labeled ApprovedFits. For more information, consult the header and comments inside `ApproveMCMCResults.m`.

# Bayes Factor Analysis
Work in progress.

# Example using test data
The repository contains a test dataset to familiarize yourself with the functionality of the MCMC fitting procedure, located in `TestScripts\TestData.mat`. This is a MATLAB structure array containing the MS2 and PP7 fluorescence signals from the 299 cells analyzed in the Liu et. al. manuscript, formatted in the fashion needed for running `TranscriptionCycleMCMC.m`.
