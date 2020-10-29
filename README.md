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
For more details on the MCMC approach, refer to Section S3 of Liu et. al. (2020). To run an MCMC fit, TranscriptionCycleMCMC.m must be executed on a dataset, which must be saved as a MATLAB .mat file. The required data format is as follows:

- there must be a structure array named <code>data</code>

Data format:

Input descriptions:
%   'fileDir', fileDir: directory of dataset files (default = root
%   directory)
%   'saveLoc', saveLoc: directory of saved MCMC results (default = root
%   directory)
%   'numParPools', numParPools: number of parallel workers to use (default
%   = 0)
%   'n_burn', n_burn: number of burn-in steps in MCMC algorithm (default = 10000)
%   'n_steps', n_steps: number of steps in MCMC algorithm including burn-in (default =
%   20000)
%   'ratePriorWidth', ratePriorWidth: standard deviation of Gaussian prior
%   for rate fluctuation term dR(t) (default = 50, see Liu et al, Section S3.1) 
%   't_start', t_start: time value to start fit at (default = 0)
%   't_end', t_end: time value to end fit at (default = Inf)
%   'loadPrevious', true/false: option to load previous inference results
%   to retain inferred elongation rate for hierarchical fit (see Liu et al,
%   Section S3.2
%   'construct', construct: option to specify a custom reporter gene that
%   must be defined in the subfunction GetFluorFromPolPos (default uses the
%   P2P-MS2-lacZ-PP7 construct from Liu et al)


Customization:


# How to define custom reporter constructs

# Approving results
# Example using test data
Description of test data
Selection of inputs
