function SS = SumofSquaresConstant_MCMCstat_FreeScaling(construct,datatype,data,x,varargin)
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
% datatype: string.
%   String specifying type of data (1 or 2-color).
% data: structure
%   Structure containing time data and fluorescence data.
% x: array
%   Array of parameters (v,R)
% varargin: input options (e.g. using dwell time or basal levels)
%   'dwell': Dwell model
%   'termination': Termination model
%   'basal': Set basal fluorescence levels
%   'PolyRate': Polynomial rate approximation of the form a_0 + a_1*t + ... +
%               a_n*t^n
%   'MeanRate': Mean rate fit with fluctuations
%
% Return
% ------
% SS: sum-of-squares residual

% Note that this assumes our sampling assumes Gaussian error with fixed
% variance/standard deviation for each measurement.

%% Extract fit parameters and data

if ~isempty(varargin)
    varargin = varargin{1};
end

R_ind_shift = 0; %Default shift in starting index for rate variable
polyrate = false; %Default is to not use the polynomial rate approximation
meanrate = false; %Default is to not use the mean rate approximation with fluctuations

if any(strcmpi(varargin,'dwell'))
    R_ind_shift = R_ind_shift + 1;
end
if any(strcmpi(varargin,'RNAPterm'))
    R_ind_shift = R_ind_shift + 1;
end
if any(strcmpi(varargin,'termination'))
    R_ind_shift = R_ind_shift + 1;
end
if any(strcmpi(varargin,'basal'))
    R_ind_shift = R_ind_shift + 2;
end
if any(strcmpi(varargin,'PolyRate'))
    polyrate = true;
end
if any(strcmpi(varargin,'MeanRate'))
    meanrate = true;
end

R_ind_start = 4 + R_ind_shift; %Starting index for rate variable

v = x(1);
ton = x(2);
A = x(3);

t = data.xdata;
fluorExp = data.ydata;

%Rate fitting approximations
if polyrate
    R = PolynomialRate(x(R_ind_start:end),t);
elseif meanrate
    R = MeanRate(x(R_ind_start:end));
else
    R = x(R_ind_start:end);
end

%% Set up Posterior function
if strcmp(datatype,'2-color')
    twocolor = 1;
else
    twocolor = 0;
end

%Calculate simulated fluorescence.
PolPos = ConstantElongationSim(v,ton,R,t);
[MS2,PP7] = GetFluorFromPolPos(construct,PolPos,varargin);
MS2 = A * MS2;

if twocolor
    fluorSim = [MS2,PP7];
else
    fluorSim = MS2;
end


%% Compute sum-of-squares
% Calculate the residuals of experimental and theoretical predictions
residuals = fluorExp - fluorSim;

% Compute the sum-of-squares function
SS = nansum(residuals.^2);
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