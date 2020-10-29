# TranscriptionCycleInference
Repository for transcription cycle MCMC inference code from Liu et. al. (2020) (https://www.biorxiv.org/content/10.1101/2020.08.29.273474v3).


Outline:
- Overview of code
- MCMCanalysis settings:
  - n_burn
  - n_steps
  - ratePriorWidth
  - t_start
  - t_end
  - loadPrevious
  - numParPool
  - construct
- GetFluorFromPolPos:
  - custom constructs
- ApproveMCMCResults
  - what to look for
- test dataset and scripts
  
ToDo:
- user option to customize initial guess and MCMC bounds/priors
- user option to load previous datasets

# Introduction
This repository contains the MATLAB code necessary to run the Markov Chain Monte Carlo (MCMC) fitting on dual-color fluorescent reporter data as described in Liu et. al. (2020). First, this code makes use of the MCMCStat package developed by Haario. et. al., found here: https://github.com/mjlaine/mcmcstat.

# Organization
The necessary scripts to run can be found in the \src folder. The primary script is TranscriptionCycleMCMC.m, which contains the master code that handles the MCMC fits. GetFluorFromPolPos.m contains the information on the reporter constructs needed to calculate fluorescent signals from the transcription cycle model, and should be modified to incorporate the reporter construct used by the user. ApproveMCMCResults.m is a manual curation script that can be used to visualized fit results and approve/reject them. The \src\dependencies folder contains sub-functions that are necessary, but should not be altered by the user.

# Guide to MCMC fitting
For more details on the MCMC approach, refer to Section S3 of Liu et. al. (2020). To run an MCMC fit, TranscriptionCycleMCMC.m must be executed on a dataset, which must be saved as a MATLAB .mat file. The required data format for each dataset is as follows: there must be a structure array named <code>data</code>. Each index of <code>data</code> corresponds to a single cell of fluorescent data. Each index (i.e. single cell) of <code>data</code> must have four fields:

- <code>time</code>: a 1 x N vector of time points of experimental acquisition
- <code>MS2</code>: a 1 x N vector of MS2 fluorescent signals
- <code>PP7</code>: a 1 x N vector of PP7 fluorescent signals
- <code>name</code>: a string used to label this dataset
  
TranscriptionCycleMCMC.m loads one or more datasets with this format and runs the MCMC inference procedure on them. To customize the fitting, the below variable arguments can be passed in:

- 'fileDir', string: directory of where to search for dataset files (default is current working directory)
- 'saveLoc', string: directory of where to save MCMC fit results (default is current working directory0
- 'numParPools', int: number of parallel pool workers to use (default = 0)
- 'n_burn', int: number of burn-in steps for MCMC algorithm (default = 10000)
- 'n_steps', int: number of steps in MCMC algorithm including burn-in (Default = 200000)
- 'ratePriorWidth', float: standard deviation of Gaussian prior for rate fluctuation term dR(t) (default = 50, see Section S3.1)
- 't_start', float: time value to start fit at (default = 0)
- 't_end', float: time value to end fit at (default = 0)
- 'loadPrevious', true/false: option to load previous inference results to retain inferred elongation rate for hierarchical fit (see Section S3.2, default = false)
- 'construct', string: option to specify a custom reporter gene that must be defined in the script GetFluorFromPolPos.m (default uses the P2P-MS2-lacZ-PP7 construct from Liu et. al. 2020)

# How to define custom reporter constructs
To use a different reporter gene than the P2P-MS2-lacZ-PP7 construct in Liu et. al. 2020, simply define a construct in the header of GetFluorFromPolPos.m. There is a template provided in the comments, where you can edit the positions of the MS2 and PP7 stem loop sequences as well as the number of stem loops per sequence. Make sure to define your reporter gene with a label, which you then pass as an argument into TranscriptionCycleMCMC.m with the variable argument 'construct'.

# Approving results
# Example using test data
Description of test data
Selection of inputs
