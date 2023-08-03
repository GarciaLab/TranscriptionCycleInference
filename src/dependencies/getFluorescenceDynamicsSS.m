function [SS,simMS2,simPP7] = getFluorescenceDynamicsSS(data,x,ElongationSegments,velocity_names,stemloops,x_stall)
%getFluorescenceDynamicsSS computes either the sum-off-squares function
%(given data and parameters, if mode==true) or simulated fluorescence traces (given
%parameters, if mode=false).

%% 1) Extract fit parameters and data

%Create interpolated time for model fitting
t = data.xdata{1}; %Time vector
t_interp = data.xdata{2}; %Interpolated time vector with even time resolution
dt = t_interp(2) - t_interp(1); %Distance between interpolated times

%Extract fluorescence values
fluorExp = data.ydata;

%Assign velocity parameters to segments
velocity_assigment = ElongationSegments.velocities;
v = zeros(1,length(velocity_assigment));
for idx = 1:length(velocity_names)
    indices = cell2mat(cellfun(@(x)(strcmpi(x,velocity_names{idx})),velocity_assigment,'UniformOutput',false));
    v(indices) = x(idx);
end

%Assign other parameters
N_vel = length(velocity_names); %Get number of velocity parameters
tau = x(N_vel+1); %Cleavage time
ton = x(N_vel+2); %Onset time
MS2_basal = x(N_vel+3); %Basal MS2 fluorescence
PP7_basal = x(N_vel+4); %Basal PP7 fluorescence
A = x(N_vel+5); %Scaling factor between MS2/PP7

if isempty(x_stall)
    R = x(N_vel+6)+ x((N_vel+7):end); %Overall initiation rate = Mean initiation rate + vector of fluctiations in initiation rate
else
    R = x(N_vel+6)+ x((N_vel+9):end); %Overall initiation rate = Mean initiation rate + vector of fluctiations in initiation rate
    ProbPremTerm = x(N_vel+7);
    tauPremTerm = x(N_vel+8);
end


%% 2) Simulate fluorescence dynamics

%If there are any negative rates, make them zero.
R(R<0) = 0;

%First index of t for which t>=ton is true.
idx_on = find(t_interp >= ton,1);

%Translate initiation rate to number of initiations
N_R = length(R);
R = R * triu(ones(N_R,N_R));
Init = floor(R*dt);
Init = Init(1:end) - [0,Init(1:end-1)];
Init(1:(idx_on-1)) = 0;
N_init = length(Init); %Note that N_init == N_R == length(t_interp)


if isempty(x_stall) % Without premature termination
    %Get possible positions of single Polymerase at equidistant time points
    %after its initiation (and effective construct length L)
    [pos,N_pos] = generatePos(dt,ElongationSegments.segments,v,tau);

    % Get fluorescence intensity of a single Polymerase at the positions pos
    [singlePolMS2,singlePolPP7,MS2_0,PP7_0] = getSinglePolFluor(pos,stemloops);

    % Get kymograph of Polymerases
    kymograph = getKymographSS(Init,N_init,N_pos,idx_on);

    % Get fluorescence dynamics
    simMS2 = singlePolMS2 * kymograph + MS2_0;
    simPP7 = singlePolPP7 * kymograph + PP7_0;
else % With deterministic premature termination model

    %1. Get possible positions of single Polymerase at equidistant time points
    %after its initiation (and effective construct length L)
    %2. Get indices of stalling times at stalling site
    [pos,N_pos,idx_stall] = generatePosStall(dt,ElongationSegments.segments,v,tau,x_stall,tauPremTerm,N_init);

    % Get fluorescence intensity of a single Polymerase at the positions pos
    [singlePolMS2,singlePolPP7,singlePolMS2Stall,singlePolPP7Stall,N_stall,MS2_0,PP7_0] = getSinglePolFluorStall(pos,stemloops,x_stall);

    % Get kymograph of Polymerases
    [kymograph,kymographPremTerm] = getKymographPremTerm(Init,N_init,N_pos,idx_on,idx_stall,N_stall,ProbPremTerm);

    % Get fluorescence dynamics
    simMS2 = singlePolMS2 * kymograph + singlePolMS2Stall * kymographPremTerm + MS2_0;
    simPP7 = singlePolPP7 * kymograph + singlePolPP7Stall * kymographPremTerm + PP7_0;
end

% Basal activity: Set fluorescence levels below detection threshold to basal level
simMS2(simMS2<MS2_basal) = MS2_basal; %MG: TODO -> Consider if it is better to first rescale and then truncate.
simPP7(simPP7<PP7_basal) = PP7_basal;

simMS2 = A * simMS2; %Rescale MS2 from PP7-based units to MS2 associated units

%Interpolate the simulated fluorescence signals back to experimental time
%resolution
simMS2 = interp1(t_interp,simMS2,t);
simPP7 = interp1(t_interp,simPP7,t);

%% 3) Compute sum-of-squares
% Calculate the scaled residuals of experimental data and theoretical predictions
residuals = ((fluorExp - [simMS2,simPP7]).^2)./fluorExp;


% Compute the sum-of-squares function
SS = nansum(residuals);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunctions
% Some functions are implemented as subfunctions to increase efficiency.
% All subfunctions ending with '*Stall' or '*PremTerm' are modifications of their
% counterparts, that do not account for premature termination.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos,N_pos] = generatePos(dt,segments,velocities_params,tau)
%generatePos evaluates the positions of a single polymerase after multiple
%equidistant time steps. It evaluates these positions by defining a
%piecwise function position(time) and assigns pos(i)=position(time(i)),
%where time(1)=0. Only those polymerases contribute to the fluorescence
%signal that have not terminated yet, implying that all positions beyond
%the effective construct length can be neglected.

L = segments(end) + velocities_params(end)*tau; %Effective length of construct (including cleavage time)
segments(end) = L; %Change end of last segment to effective length
N_seg = length(segments); %Number of segments

% Get intervals for piecewise definition of the polymerase position as a
% function of time
t_seg = zeros(1,1+N_seg); %Initialize and add time=0
t_seg(1) = 0;
t_seg(2) = segments(1)/velocities_params(1);
if N_seg>1
    for k=3:(1+N_seg)
        t_seg(k) = t_seg(k-1) + (segments(k-1)-segments(k-2))/velocities_params(k-1);
    end
end

times = 0:dt:t_seg(end); %Get interpolated acquisition times, assuming single Polymerase initiation at t=0
N_pos = length(times); %Length of the position vector
pos = zeros(1,N_pos); %Initialize positions of the Polymerase for all interpolated acquisition times

% Assign position values of a single Polymerase at equidistant times after initiation
idx = (times<=t_seg(2));
pos(idx) = velocities_params(1)*times(idx);
if N_seg>1
    for k=3:(1+N_seg)
        idx = and(times<=t_seg(k),times>t_seg(k-1));
        pos(idx) = segments(k-2) + velocities_params(k-1)*(times(idx)-t_seg(k-1));
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [singlePolMS2,singlePolPP7,MS2_0,PP7_0] = getSinglePolFluor(pos,stemloops)
%Get a row-vectors of MS2 and PP7 fluorescence intensity, corresponding to a
%single Polymerase at positions pos. These vectors can be multiplied
%(matrix product) with the rows of the kymograph matrix (space-time plot of
%Polymerase numbers) to get thefluorescence intensities at all times.

%Number of positions
N_pos = length(pos);

%Initialize
MS2_start = stemloops.MS2_start;
PP7_start = stemloops.PP7_start;
MS2_end = stemloops.MS2_end;
PP7_end = stemloops.PP7_end;
MS2_loopn = stemloops.MS2_loopn;
PP7_loopn = stemloops.PP7_loopn;
MS2_0 = stemloops.MS2_0;
PP7_0 = stemloops.PP7_0;
singlePolMS2 = zeros(1,N_pos);
singlePolPP7 = zeros(1,N_pos);

% Total number of loops (+ constitutively bound fluorophores,i.e., bound to the DNA not RNA)
MS2_Sum = sum(MS2_loopn) + MS2_0;
PP7_Sum = sum(PP7_loopn) + PP7_0;

% Rescaling of number of constitutively bound fluorophores to intensity units
MS2_0 = MS2_0/MS2_Sum; 
PP7_0 = PP7_0/PP7_Sum;

%Initialize without constitutively bound fluorophores
MS2 =  0;
PP7 = 0;
singlePolMS2(pos<=MS2_start(1)) = MS2;
singlePolPP7(pos<=PP7_start(1)) = PP7;

% Loop over MS2 subsequences (except last one)
if length(MS2_start)>1
    for i = 1:(length(MS2_start)-1)
        idx = and(pos>MS2_start(i),pos<=MS2_end(i));
        MS2_Fluorval = MS2_loopn(i)/MS2_Sum; % Fraction of MS2 loops in i-th subsequence
        singlePolMS2(idx) = MS2 + (pos(idx)- MS2_start(i)).*MS2_Fluorval./(MS2_end(i) - MS2_start(i));
        MS2 = MS2 + MS2_Fluorval;
        idx = and(pos>MS2_end(i),pos<=MS2_start(i+1));
        singlePolMS2(idx) = MS2;
    end
end

%Catch up on last subsequence
idx = and(pos>MS2_start(end),pos<=MS2_end(end));
MS2_Fluorval = MS2_loopn(end)/MS2_Sum; % Fraction of MS2 loops in i-th subsequence
singlePolMS2(idx) = MS2 + (pos(idx)- MS2_start(end)).*MS2_Fluorval./(MS2_end(end) - MS2_start(end));
MS2 = MS2 + MS2_Fluorval;
singlePolMS2(pos>MS2_end(end)) = MS2;


% Loop over PP7 subsequences
if length(PP7_start)>1
    for i = 1:(length(PP7_start)-1)
        idx = and(pos>PP7_start(i),pos<=PP7_end(i));
        PP7_Fluorval = PP7_loopn(i)/PP7_Sum; % Fraction of MS2 loops in i-th subsequence
        singlePolPP7(idx) = PP7 + (pos(idx)- PP7_start(i)).*PP7_Fluorval./(PP7_end(i) - PP7_start(i));
        PP7 = PP7 + PP7_Fluorval;
        idx = and(pos>PP7_end(i),pos<=PP7_start(i+1));
        singlePolPP7(idx) = PP7;
    end
end

%Catch up on last subsequence
idx = and(pos>PP7_start(end),pos<=PP7_end(end));
PP7_Fluorval = PP7_loopn(end)/PP7_Sum; % Fraction of MS2 loops in i-th subsequence
singlePolPP7(idx) = PP7 + (pos(idx)- PP7_start(end)).*PP7_Fluorval./(PP7_end(end) - PP7_start(end));
PP7 = PP7 + PP7_Fluorval;
singlePolPP7(pos>PP7_end(end)) = PP7;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kymograph = getKymographSS(Init,N_t,N_pos,idx_on)
%Converts vector of the number of initiating polymerases at different
%interpolated times into a discretized kymograph (position x time)

kymograph = zeros(N_pos,N_t); %Initialize kymograph (positions x times)

idx_1 = (idx_on-1)*N_pos +1; %Get index of first non-zero element of the kymograph

if (N_t-idx_on+1-N_pos)>=0 %Account for the possibility of late t_on or short observation times
    %Start with case of early t_on or long observation times
    idx = idx_1:(N_pos+1):(N_pos+idx_on-1)*N_pos; %Get indices of the (off-)diagonal starting at element idx_1

    % Loop over all (off-)diagonals of full length
    for k=0:(N_t-idx_on+1-N_pos)
        kymograph(idx) = Init(idx_on+k); %Set (off-)diagonals of full length to constant value
        idx=idx+N_pos; %Get indices of next off-diagonal
    end

    % Loop over all (off-)diagonals of reduced length
    for k=1:N_pos-1
        kymograph(idx(1:(end-k))) = Init(N_t-N_pos+1+k); %Set (off-)diagonals of reduced length to a constant value
        idx=idx+N_pos; %Get indices of next off-diagonal
    end
else
    %For some chain elements the observation time may be shorter than the time the first poymerase needs to traverse the construct.
    idx = idx_1:(N_pos+1):(N_t*N_pos);
    idx = idx(1:(N_t-idx_on+1)); %Truncate index set of first (off-)diagonal to the correct length

    % Loop over all (off-)diagonals of reduced length
    for k=0:(N_t-idx_on)
        kymograph(idx(1:(end-k))) = Init(idx_on+k); %Set (off-)diagonals of reduced length to a constant value
        idx=idx+N_pos; %Get indices of next off-diagonal
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos,N_pos,idx_stall] = generatePosStall(dt,segments,velocities_params,tau,x_stall,tauPremTerm,N_t)
%generatePosStall evaluates the positions of a single polymerase after multiple
%equidistant time steps. It evaluates these positions by defining a
%piecwise function position(time) and assigns pos(i)=position(time(i)),
%where time(1)=0. Only those polymerases contribute to the fluorescence
%signal that have not terminated yet, implying that all positions beyond
%the effective construct length can be neglected.
%
%Additionally to the function generatePos it returns the array of indices
%idx_stall, which contains in each row the lower and upper time indices, between which a
%polymerase, which is about to terminate prematurely, stalls.

L = segments(end) + velocities_params(end)*tau; %Effective length of construct (including cleavage time)
segments(end) = L; %Change end of last segment to effective length
N_seg = length(segments); %Number of segments

% Get intervals for piecewise definition of the polymerase position as a
% function of time
t_seg = zeros(1,1+N_seg); %Initialize and add time=0
t_seg(1) = 0;
t_seg(2) = segments(1)/velocities_params(1);
if N_seg>1
    for k=3:(1+N_seg)
        t_seg(k) = t_seg(k-1) + (segments(k-1)-segments(k-2))/velocities_params(k-1);
    end
end

times = 0:dt:t_seg(end); %Get interpolated acquisition times, assuming single Polymerase initiation at t=0
N_pos = length(times); %Length of the position vector
pos = zeros(1,N_pos); %Initialize positions of the Polymerase for all interpolated acquisition times

% Assign position values of a single Polymerase at equidistant times after initiation
idx = (times<=t_seg(2));
pos(idx) = velocities_params(1)*times(idx);
if N_seg>1
    for k=3:(1+N_seg)
        idx = and(times<=t_seg(k),times>t_seg(k-1));
        pos(idx) = segments(k-2) + velocities_params(k-1)*(times(idx)-t_seg(k-1));
    end
end

% Evaluating the inverse of the function position(time) at x_stall(l) to get
% the times when the polymerase arrives at x_stall(l). Then get extract the
% time indices between which the polymerase stalls before terminating prematurely.
idx_stall = zeros(length(x_stall),2);
for l=1:length(x_stall)
    idx_seg = find(segments>x_stall(l),1);
    if idx_seg==1
        t_stall = x_stall(l)/velocities_params(1); %Evaluate in the case of just a single velocity parameter
    else
        t_stall = (x_stall(l) - segments(idx_seg-1))/velocities_params(idx_seg) + t_seg(idx_seg) ;
    end
    idx_stall(l,2) = min(N_t,find((t_stall + tauPremTerm)>=times,1,'last')); %Get last index satisfying the condition (last time before premature termination)
    idx_stall(l,1) = min(N_t,find(t_stall<=times,1)); %Get first index satisfying the condition (last stalling time)
    %The min is necessary for idx_stall(l,2), since the number of positions
    %can be larger than the number of interpolated times.
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [singlePolMS2,singlePolPP7,singlePolMS2_Stall,singlePolPP7_Stall,N_stall,MS2_0,PP7_0] = getSinglePolFluorStall(pos,stemloops,x_stall)
%Get a row-vectors of MS2 and PP7 fluorescence intensity, corresponding to a
%single Polymerase at positions pos and positions x_stall. These vectors can be multiplied
%(matrix product) with the rows of a kymograph matrix (space-time plot of
%Polymerase numbers) to get thefluorescence intensities at all times.

%Number of positions
N_pos = length(pos);
N_stall = length(x_stall);

%Initialize
MS2_start = stemloops.MS2_start;
PP7_start = stemloops.PP7_start;
MS2_end = stemloops.MS2_end;
PP7_end = stemloops.PP7_end;
MS2_loopn = stemloops.MS2_loopn;
PP7_loopn = stemloops.PP7_loopn;
MS2_0 = stemloops.MS2_0;
PP7_0 = stemloops.PP7_0;
singlePolMS2 = zeros(1,N_pos);
singlePolPP7 = zeros(1,N_pos);
singlePolMS2_Stall = zeros(1,N_stall);
singlePolPP7_Stall = zeros(1,N_stall);

% Total number of loops (+ constitutively bound fluorophores,i.e., bound to the DNA not RNA)
MS2_Sum = sum(MS2_loopn) + MS2_0;
PP7_Sum = sum(PP7_loopn) + PP7_0;

% Rescaling of number of constitutively bound fluorophores to intensity units
MS2_0 = MS2_0/MS2_Sum; 
PP7_0 = PP7_0/PP7_Sum;

%Initialize without constitutively bound fluorophores
MS2 =  0;
PP7 = 0;
singlePolMS2(pos<=MS2_start(1)) = MS2;
singlePolPP7(pos<=PP7_start(1)) = PP7;
singlePolMS2_Stall(x_stall<=MS2_start(1)) = MS2;
singlePolPP7_Stall(x_stall<=PP7_start(1)) = PP7;

% Loop over MS2 subsequences (except last one)
if length(MS2_start)>1
    for i = 1:(length(MS2_start)-1)
        idx = and(pos>MS2_start(i),pos<=MS2_end(i));
        idx_stall = and(x_stall>MS2_start(i),x_stall<=MS2_end(i));
        MS2_Fluorval = MS2_loopn(i)/MS2_Sum; % Fraction of MS2 loops in i-th subsequence
        singlePolMS2(idx) = MS2 + (pos(idx)- MS2_start(i)).*MS2_Fluorval./(MS2_end(i) - MS2_start(i));
        singlePolMS2_Stall(idx_stall) = MS2 + (x_stall(idx_stall)- MS2_start(i)).*MS2_Fluorval./(MS2_end(i) - MS2_start(i));
        MS2 = MS2 + MS2_Fluorval;
        idx = and(pos>MS2_end(i),pos<=MS2_start(i+1));
        idx_stall = and(x_stall>MS2_end(i),x_stall<=MS2_start(i+1));
        singlePolMS2(idx) = MS2;
        singlePolMS2_Stall(idx_stall) = MS2;
    end
end

%Catch up on last subsequence
idx = and(pos>MS2_start(end),pos<=MS2_end(end));
idx_stall = and(x_stall>MS2_start(end),x_stall<=MS2_end(end));
MS2_Fluorval = MS2_loopn(end)/MS2_Sum; % Fraction of MS2 loops in i-th subsequence
singlePolMS2(idx) = MS2 + (pos(idx)- MS2_start(end)).*MS2_Fluorval./(MS2_end(end) - MS2_start(end));
singlePolMS2_Stall(idx_stall) = MS2 + (x_stall(idx_stall)- MS2_start(end)).*MS2_Fluorval./(MS2_end(end) - MS2_start(end));
MS2 = MS2 + MS2_Fluorval;
singlePolMS2(pos>MS2_end(end)) = MS2;
singlePolMS2_Stall(x_stall>MS2_end(end)) = MS2;


% Loop over PP7 subsequences
if length(PP7_start)>1
    for i = 1:(length(PP7_start)-1)
        idx = and(pos>PP7_start(i),pos<=PP7_end(i));
        idx_stall = and(x_stall>PP7_start(i),x_stall<=PP7_end(i));
        PP7_Fluorval = PP7_loopn(i)/PP7_Sum; % Fraction of PP7 loops in i-th subsequence
        singlePolPP7(idx) = PP7 + (pos(idx)- PP7_start(i)).*PP7_Fluorval./(PP7_end(i) - PP7_start(i));
        singlePolPP7_Stall(idx_stall) = PP7 + (x_stall(idx_stall)- PP7_start(i)).*PP7_Fluorval./(PP7_end(i) - PP7_start(i));
        PP7 = PP7 + PP7_Fluorval;
        idx = and(pos>PP7_end(i),pos<=PP7_start(i+1));
        idx_stall = and(x_stall>PP7_end(i),x_stall<=PP7_start(i+1));
        singlePolPP7(idx) = PP7;
        singlePolPP7_Stall(idx_stall) = PP7;
    end
end

%Catch up on last subsequence
idx = and(pos>PP7_start(end),pos<=PP7_end(end));
idx_stall = and(x_stall>PP7_start(end),x_stall<=PP7_end(end));
PP7_Fluorval = PP7_loopn(end)/PP7_Sum; % Fraction of PP7 loops in i-th subsequence
singlePolPP7(idx) = PP7 + (pos(idx)- PP7_start(end)).*PP7_Fluorval./(PP7_end(end) - PP7_start(end));
singlePolPP7_Stall(idx_stall) = PP7 + (x_stall(idx_stall)- PP7_start(end)).*PP7_Fluorval./(PP7_end(end) - PP7_start(end));
PP7 = PP7 + PP7_Fluorval;
singlePolPP7(pos>PP7_end(end)) = PP7;
singlePolPP7_Stall(x_stall>PP7_end(end)) = PP7;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kymograph,kymographPremTerm] = getKymographPremTerm(Init,N_t,N_pos,idx_on,idx_stall,N_stall,ProbPremTerm)
%Converts vector of the number of initiating polymerases at different
%interpolated times and the specifications of premature termination into two discretized
%kymographs (position x time), one kymograph for all progressive
%polymerases and one kymograph for all stalling polymerases at stalling
%sites.
%
% Premature termination dynamics:
% At each stalling site the number of polymerases (in units of PP7) is
% deterministically reduced by the constant factor ProbPremTerm.
% kymographPremTerm is the "kymograph" restricted to the stalling sites, i.e.,
% it has dimensions (N_stall x N_t).

SurviveProportion = ((1-ProbPremTerm).^(1:N_stall))'; % Generate vector of proportions of Polymerases that survive after l sites
SurvivorVecs = repmat(Init,N_stall,1).*SurviveProportion; % Generate stack of vectors of the surving polymerases after the individual stalling sites
PremTermVecs = repmat(Init,N_stall,1)-SurvivorVecs; % Get analogous stack of vectors for polymerases that stall and terminate prematurely


kymograph = zeros(N_pos,N_t); %Initialize kymograph (positions x times)
kymographPremTerm = zeros(N_stall,N_t); %Initialize kymographPremTerm (stalling positions x times)
PremTermVec2Kymograph = zeros(N_t,N_t,N_stall); %Initialize conversion matrices from PremTermVecs to rows of kymographPremTerm

idx_1 = (idx_on-1)*N_pos +1; %Get index of first non-zero element of the kymograph

if (N_t-idx_on+1-N_pos)>=0 %Account for the possibility of late t_on or short observation times
    %Start with case of early t_on or long observation times
    idx = idx_1:(N_pos+1):(N_pos+idx_on-1)*N_pos; %Get indices of the (off-)diagonal starting at element idx_1

    % Loop over all (off-)diagonals of full length
    for k=0:(N_t-idx_on+1-N_pos)
        kymograph(idx(1:idx_stall(1,1)-1)) = Init(idx_on+k); %Set (off-)diagonals up to first stalling site to constant value
        for l=1:N_stall
            kymograph(idx(idx_stall(l,1):end)) = SurvivorVecs(l,idx_on+k); %Set (off-)diagonals after l-th stalling site to constant value
        end
        idx=idx+N_pos; %Get indices of next off-diagonal
    end

    % Loop over all (off-)diagonals of reduced length
    for k=1:N_pos-1
        varidx = idx(1:(end-k)); %Define truncated index set
        kymograph(varidx(1:min(end,idx_stall(1,1)-1))) = Init(N_t-N_pos+1+k); %Set (off-)diagonals up to first stalling site to constant value
        for l=1:N_stall
            kymograph(varidx(idx_stall(l,1):end)) = SurvivorVecs(l,N_t-N_pos+1+k); %Set (off-)diagonals after l-th stalling site to constant value
        end
        idx=idx+N_pos; %Get indices of next off-diagonal
    end
else
    %For some chain elements the observation time may be shorter than the time the first poymerase needs to traverse the construct.
    idx = idx_1:(N_pos+1):(N_t*N_pos);
    idx = idx(1:(N_t-idx_on+1)); %Truncate index set of first (off-)diagonal to the correct length

    % Loop over all (off-)diagonals of reduced length
    for m=0:(N_t-idx_on)
        varidx = idx(1:(end-m)); %Define truncated index set
        kymograph(varidx(1:min(end,idx_stall(1,1)-1))) = Init(idx_on+m); %Set (off-)diagonals up to first stalling site to constant value
        for l=1:N_stall
            kymograph(varidx(idx_stall(l,1):end)) = SurvivorVecs(l,idx_on+m); %Set (off-)diagonals after l-th stalling site to constant value
        end
        idx=idx+N_pos; %Get indices of next off-diagonal
    end
end


% Calculate kymographPremTerm via matrix multiplication for each site l
for l=1:N_stall
    %Get conversion matrix for l-th stalling site
    for m=idx_stall(l,1):idx_stall(l,2)
        PremTermVec2Kymograph(:,:,l) = PremTermVec2Kymograph(:,:,l) + diag(ones(N_t-m+1,1),m-1); %Add diagonals to account for polymerases that have reached the stalling site and stall there for the time tauPremTerm
    end
    PremTermVec2Kymograph(1:(idx_on-1),:,l) = 0; % Remove the contribution of initiations before t_on

    % Calculate time series of polymerase number at l-th stalling site
    kymographPremTerm(l,:) = PremTermVecs(l,:) * PremTermVec2Kymograph(:,:,l);
end

end