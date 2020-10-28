function SS = SumofSquaresFunction_TranscriptionCycleMCMC(construct,data,x)
% Log likelihood probability function for constant model elongation rate MCMC fitting.
% Currently the model infers values for elongation rate v, noise term D,
% loading rate R, and measurement noise sigma. This is for use in the
% MCMCstat package and returns the sum-of-squares function. We assume that
% the calibration factor between MS2 and PP7 is unknown and use it as a
% free parameter.

% Parameters
% ----------
% construct: string.
%   String specifying the construct used.
% data: structure
%   Structure containing time data and fluorescence data.
% x: array
%   Vector of fit parameters (v,tau,t_on,MS2_basal,PP7_basal,A,R,dR)

% Return
% ------
% SS: sum-of-squares residual

% Note that this assumes our sampling assumes Gaussian error with fixed
% variance/standard deviation for each measurement.

%% Extract fit parameters and data

%Create interpolated time for model fitting
t = data.xdata; %Time vector
dt = mean(t(2:end)-t(1:(end-1))); %Average time resolution of the dataset for model usage
t_interp = t(1):dt:t(end); %Interpolated time vector with even time resolution

%Extract fluorescence values
fluorExp = data.ydata;

v = x(1); %Elongation rate
tau = x(2); %Cleavage rate
ton = x(3); %Onset time
MS2_basal = x(4); %Basal MS2 fluorescence
PP7_basal = x(5); %Basal PP7 fluorescence
A = x(6); %Scaling factor between MS2/PP7
R = x(7); %Mean initiation rate
dR = x(8:end); %Fluctiations in initiation rate

%Calulate overall initiation rate
R_full = R + dR;

%% Set up posterior function for MCMC
%Calculate simulated fluorescence.
PolPos = ConstantElongationSim(v,ton,R_full,t_interp); %Simulated polymerase positions
[MS2,PP7] = GetFluorFromPolPos(construct,PolPos,v,tau,MS2_basal,PP7_basal); %Simulated fluorescences
MS2 = A * MS2; %MS2 after scaling factor

%Interpolate the simulated fluorescence signals back to experimental time
%resolution
MS2 = interp1(t_interp,MS2,t);
PP7 = interp1(t_interp,PP7,t);
fluorSim = [MS2,PP7]; %Simulated fluorescences

%% Compute sum-of-squares
% Calculate the residuals of experimental and theoretical predictions
residuals = fluorExp - fluorSim;

% Compute the sum-of-squares function
SS = nansum(residuals.^2);
end