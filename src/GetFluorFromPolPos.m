function [MS2,PP7] = GetFluorFromPolPos(construct,PolPos,v,tau,MS2_basal,PP7_basal)
%Calculates an MS2 and PP7 signal given a simulated position matrix of Pol
%II molecules over time on a gene. If a Pol II is past the end of a loop
%sequence, give it a unit fluorescence (or # count of Pol II). If it is
%partially done through a loop sequence, give it a fractional unit of
%fluorescence.

%Inputs:
%   construct: string defining the construct used
%   PolPos: matrix giving positions of Pol II molecules over time (time x
%   space)
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
%    MS2_start = <POSITION OF START OF MS2 LOOPS>; %Position of start of MS2 loops
%    MS2_end = <POSITION OF END OF MS2 LOOPS>; %Position of end of MS2 loops
%    PP7_start = <POSITION OF START OF PP7 LOOPS>; %Position of start of PP7 loops
%    PP7_end = <POSITION OF END OF PP7 LOOPS>; %Position of end of PP7 loops
%    
%    MS2_loopn = <NUMBER OF MS2 LOOPS>; %Number of MS2 loops
%    PP7_loopn = <NUMBER OF PP7 LOOPS>; %Number of PP7 loops
%    
%    MS2 = 0;
%    PP7 = 0;
end
%Calculate MS2 signal by making a fluorescence map
for i = 1:length(MS2_start)
    MS2_fluorval = MS2_loopn(i)/24;
    MS2map = zeros(size(PolPos));
    MS2map(and(PolPos > MS2_end(i), PolPos <L_MS2)) = MS2_fluorval; %Pol II's that are past the MS2 sequence.
    frac_ind = and(PolPos > MS2_start(i), PolPos < MS2_end(i)); %Indices for partially transcribed loops.
    MS2map(frac_ind) = (PolPos(frac_ind) - MS2_start(i)).*MS2_fluorval./(MS2_end(i) - MS2_start(i));

    MS2 = MS2 + sum(MS2map,2)';
    
    %Basal activity
    MS2(MS2<MS2_basal) = MS2_basal;

    %Calculate PP7 signal by making a fluorescence map
    PP7_fluorval = PP7_loopn(i)/24;
    PP7map = zeros(size(PolPos));
    PP7map(and(PolPos > PP7_end(i), PolPos <L_PP7)) = PP7_fluorval; %Pol II's that are past the MS2 sequence.
    frac_ind = and(PolPos > PP7_start(i), PolPos < PP7_end(i)); %Indices for partially transcribed loops.
    PP7map(frac_ind) = (PolPos(frac_ind) - PP7_start(i)).*PP7_fluorval./(PP7_end(i) - PP7_start(i));

    PP7 = PP7 + sum(PP7map,2)';
    
    %Basal activity
    PP7(PP7<PP7_basal) = PP7_basal;
end
    
end