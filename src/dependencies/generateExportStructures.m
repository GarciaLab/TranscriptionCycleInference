function [fields_params,fields_MCMCchain,fields_MCMCresults,MCMCchain,MCMCresults,MCMCplot] = generateExportStructures(N,velocity_names,x_drop)
%generateExportStructures generates (a) lists of names and (b) structures
%for saving the results of the MCMC analysis of N individual cells.
%
%fields_params: list of all parameters, expect of 's2'
%fields_MCMCchain: cell array of field names for the export structure MCMCchain
%fields_MCMCresults: cell array of field names for the export structure MCMCresults
%MCMCchain: Empty (1xN) structure array with fields fields_MCMCchain
%MCMCresults: Empty (1xN) structure array with fields fields_MCMCresults
%MCMCplot: Empty (1xN) structure array with fields 't_plot', 'MS2_plot', 'PP7_plot', 'simMS2', 'simPP7'

%Generate list of parameters (either with or without drop-off)
if isempty(x_drop)
    params = [velocity_names,{'tau','ton','MS2_basal','PP7_basal','A','R','dR','s2'}];
else
    params = [velocity_names,{'tau','ton','MS2_basal','PP7_basal','A','R','ProbDrop','tauDrop','dR','s2'}];
    %Add additional parameters 'ProbDrop' (drop-off probability) and
    %'tauDrop' (drop-off (cleavage) time)
end

% Generate all fields
fields_MCMCchain = cellfun(@(x) strcat(x, '_chain'), params, 'UniformOutput', false); %Add suffix '_chain' to all parameter names
fields_MCMCchain{end} = 's2chain'; %Change 's2_chain' to 's2chain' (for consistency with previous version)
varparams = params; varparams{end} = 'sigma'; %Change 's2' to 'sigma' (for consistency with previous version)
fields_MCMCresults = [cellfun(@(x) strcat('mean_',x), varparams, 'UniformOutput', false),...
    cellfun(@(x) strcat('sigma_',x), varparams, 'UniformOutput', false),{'cell_index','ApprovedFits','FittedConstruct'}]; %Add prefixes 'mean_' and 'sigma_' to all parameter names. Further add the fields 'cell_index' and 'ApprovedFits'.

%Generate empty structures
MCMCchain = cell2struct(cell(1,length(fields_MCMCchain)),fields_MCMCchain,2);
MCMCresults = cell2struct(cell(1,length(fields_MCMCresults)),fields_MCMCresults,2);
MCMCplot = struct('t_plot',{},'MS2_plot',{},'PP7_plot',{},...
    'simMS2',{},'simPP7',{});

%Extend structures to arrays (important for parfor-loop)
MCMCchain(N).s2chain=[];
MCMCresults(N).cell_index=[];
MCMCplot(N).t_plot=[];

% List of parameter names without s2
fields_params = params(1:end-1);