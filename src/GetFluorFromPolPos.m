function [MS2,PP7] = GetFluorFromPolPos(construct,PolPos,v,tau,MS2_basal,PP7_basal)
%Calculates an MS2 and PP7 signal given a simulated position matrix of Pol
%II molecules over time on a gene. If a Pol II is past the end of a loop
%sequence, give it a unit fluorescence (or # count of Pol II). If it is
%partially done through a loop sequence, give it a fractional unit of
%fluorescence.

%Inputs:
%   construct: string defining the construct used
%   PolPos: matrix giving positions of Pol II molecules over time (time x
%   polymerase index -> space)
%   v: elongation rate
%   tau: cleavage time (i.e. effective additional gene length)
%   MS2_basal: basal MS2 fluorescence
%   PP7_basal: basal PP7 fluorescence

% Define construct parameters.
if strcmp(construct,'P2P-MS2v5-LacZ-PP7v4')
    L_MS2 = 6.626+tau*v; %Length in kb of construct (incorporating dwell time)
    L_PP7 = 6.626+tau*v; %Length in kb of construct (incorporating dwell time)
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 4.292; %Position of start of PP7 loops
    PP7_end = 5.758; %Position of end of PP7 loops
    
    MS2_loopn = 24; %Number of MS2 loops 
    PP7_loopn = 24; %Number of PP7 loops 
    
    MS2 = 0;
    PP7 = 0;

%Define your own custom construct by modifying the above values and
%changing the appropriate values (see template below)
%    L_MS2 = <LENGTH OF REPORTER GENE>+tau*v; %Length in kb of construct (incorporating dwell time)
%    L_PP7 = <LENGTH OF REPORTER GENE>+tau*v; %Length in kb of construct (incorporating dwell time)
%
%    % The following positions and numbers are vectors, if the construct
%    has separated stem loop sequences; the start/end position vectors need
%    to be sorted, such that the coordinate values are increasing.
%    MS2_start = <POSITION OF START OF MS2 LOOPS>; %Position of start of MS2 loops (in general a vector)
%    MS2_end = <POSITION OF END OF MS2 LOOPS>; %Position of end of MS2
%    loops (in general a vector)
%    PP7_start = <POSITION OF START OF PP7 LOOPS>; %Position of start of
%    PP7 loops (in general a vector)
%    PP7_end = <POSITION OF END OF PP7 LOOPS>; %Position of end of PP7
%    loops (in general a vector)
%    
%    MS2_loopn = <NUMBER OF MS2 LOOPS>; %Number of MS2 loops (in general a vector)
%    PP7_loopn = <NUMBER OF PP7 LOOPS>; %Number of PP7 loops (in general a vector)
%    
%    MS2 = 0;
%    PP7 = 0;
end

%% Calculate MS2/PP7 signal by making a fluorescence map
% Initialize maps
MS2map = zeros(size(PolPos));
PP7map = zeros(size(PolPos));

% Total number of loops
MS2_LoopSum = sum(MS2_loopn);
PP7_LoopSum = sum(PP7_loopn);

% Loop over MS2 subsequences
for i = 1:length(MS2_start)
    MS2_Fluorval = MS2_loopn(i)/MS2_LoopSum; % Fraction of MS2 loops in i-th subsequence
    full_ind = and(PolPos > MS2_end(i), PolPos <= L_MS2); %Indices past the i-th MS2 subsequence
    MS2map(full_ind) = MS2map(full_ind) +  MS2_Fluorval; % Add contribution of the i-th MS2 subsequence to all polymerase positions past that subsequence.
    frac_ind = and(PolPos >= MS2_start(i), PolPos <= MS2_end(i)); %Indices for partially transcribed loops of i-th subsequence.
    MS2map(frac_ind) = MS2map(frac_ind) + (PolPos(frac_ind) - MS2_start(i)).*MS2_Fluorval./(MS2_end(i) - MS2_start(i)); % Add contribution of the partially transcribed i-th MS2 subsequence to all polymerase positions within that subsequence.
end
    
% Loop over PP7 subsequences
for i = 1:length(PP7_start)
    PP7_Fluorval = PP7_loopn(i)/PP7_LoopSum; % Fraction of PP7 loops in i-th subsequence
    full_ind = and(PolPos > PP7_end(i), PolPos <= L_PP7); %Indices past the i-th PP7 subsequence
    PP7map(full_ind) = PP7map(full_ind) + PP7_Fluorval; % Add contribution of the i-th PP7 subsequence to all polymerase positions past that subsequence.
    frac_ind = and(PolPos >= PP7_start(i), PolPos <= PP7_end(i)); %Indices for partially transcribed loops within the 
    PP7map(frac_ind) = PP7map(frac_ind) + (PolPos(frac_ind) - PP7_start(i)).*PP7_Fluorval./(PP7_end(i) - PP7_start(i)); % Add contribution of the partially transcribed i-th PP7 subsequence to all polymerase positions within that subsequence.
end

% Cumulative fluorescence of all polymerases
MS2 = MS2 + sum(MS2map,2)';
PP7 = PP7 + sum(PP7map,2)';

% Basal activity: Set fluorescence levels below detection threshold to basal level
MS2(MS2<MS2_basal) = MS2_basal;
PP7(PP7<PP7_basal) = PP7_basal;

end