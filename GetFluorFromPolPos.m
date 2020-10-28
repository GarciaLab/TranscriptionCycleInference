function [MS2,PP7] = GetFluorFromPolPos(construct,PolPos,varargin)
%Calculates an MS2 and PP7 signal given a simulated position matrix of Pol
%II molecules over time on a gene. If a Pol II is past the end of a loop
%sequence, give it a unit fluorescence (or # count of Pol II). If it is
%partially done through a loop sequence, give it a fractional unit of
%fluorescence.

%Inputs:
%   construct: string defining the construct used
%   PolPos: matrix giving positions of Pol II molecules over time (time x
%   space)
%   varargin: various input options
%       'dwell': include possibility of cleavage dwell time at end of gene
%                follow this command with the dwell time and elongation
%                rate
%       'RNAPterm': include RNAP termination time (for the
%               tub3'UTR construct's 3' PP7 signal). Follow this command with the
%               RNAP termination time (assumes that the cleavage dwell time is
%               already given with the elongation rate)
%       'termination': include fixed termination rate at end of gene,
%                      incorporated as a flux loss term in fluorescence, to model the
%                      effects of finite termination rate (in contrast to e.g. a fixed
%                      dwell time per polymerase). Follow this command with the
%                      termination rate (in units of fluorescence/time), then the elongation rate.
%       'basal': include a basal fluorescence level for each signal (this
%                is the noise floor). Follow this command with the basal number of
%                RNAP molecules for the MS2, then the basal number of RNAP for the
%                PP7. For fluorescence levels below the basal level, replace the
%                signal with the basal level.

if ~isempty(varargin)
    varargin = varargin{1};
end

t_dwell = 0; %Default is to have no dwell time
t_term = 0; %Default is to have no RNAP termination time
r_termination = Inf; %Default is to have instantaneous termination rate.
v_elong = 0;
MS2_basal = 0; %Default is to have no basal activity
PP7_basal = 0;

for i=1:length(varargin)
    if strcmpi(varargin{i},'dwell')
        t_dwell = varargin{i+1};
        v_elong = varargin{i+2};
    end
    if strcmpi(varargin{i},'RNAPterm')
        t_term = varargin{i+1};
    end
    if strcmpi(varargin{i},'termination')
        r_termination = varargin{i+1};
        v_elong = varargin{i+2};
    end
    if strcmpi(varargin{i},'basal')
        MS2_basal = varargin{i+1};
        PP7_basal = varargin{i+2};
    end
end

% Define construct parameters.
if strcmp(construct,'P2P-MS2v5-LacZ-PP7v4')
    L_MS2 = 6.626+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    L_PP7 = 6.626+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 4.292; %Position of start of PP7 loops
    PP7_end = 5.758; %Position of end of PP7 loops
    
    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 24; %Number of PP7 loops
    
    MS2 = 0;
    PP7 = 0;
elseif strcmp(construct,'P2P-MS2v5-LacZ-601-PP7v4')
    L_MS2 = 6.626+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    L_PP7 = 6.626+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 4.292; %Position of start of PP7 loops
    PP7_end = 5.758; %Position of end of PP7 loops
    
    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 24; %Number of PP7 loops
    
    MS2 = 0;
    PP7 = 0;
elseif strcmp(construct,'P2P-MS2v5-LacZshort-PP7v4')
    L_MS2 = 3.906+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    L_PP7 = 3.906+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 1.566; %Position of start of PP7 loops
    PP7_end = 3.032; %Position of end of PP7 loops
    
    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 24; %Number of PP7 loops
    
    MS2 = 0;
    PP7 = 0;
elseif strcmp(construct,'P2P-MS2v5-LacZlong-PP7v4')
    L_MS2 = 9.440+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    L_PP7 = 9.440+t_dwell*v_elong; %Length in kb of construct (incorporating dwell time)
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 7.100; %Position of start of PP7 loops
    PP7_end = 8.566; %Position of end of PP7 loops
    
    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 24; %Number of PP7 loops
    
    MS2 = 0;
    PP7 = 0;
elseif strcmp(construct,'P2P-MS2v5-LacZ')
    L_MS2 = 5.160+t_dwell*v_elong; %Length in kb of construct
    L_PP7 = 5.160+t_dwell*v_elong; %Length in kb of construct
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = L_PP7; %Position of start of PP7 loops
    PP7_end = L_PP7; %Position of end of PP7 loops
    
    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 0; %Number of PP7 loops
    
    MS2 = 0;
    PP7 = 0;
elseif strcmp(construct,'P2P-MS2v5-Tub3UTR-PP7v4-Tub3UTR')
    L_MS2 = 4.21+t_dwell*v_elong; %Effective length of construct for MS2(incorporating dwell time)
    L_PP7 = 4.21+t_dwell*v_elong+t_term*v_elong; %Effective length of construct for PP7 (dwell + termination time)
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 4.217; %Position of start of PP7 loops
    PP7_end = 5.683; %Position of end of PP7 loops
    
    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 24; %Number of PP7 loops
    
    MS2 = 0;
    PP7 = 0;
end
%Calculate MS2 signal by making a fluorescence map
for i = 1:length(MS2_start)
    MS2_fluorval = MS2_loopn(i)/24;
    MS2map = zeros(size(PolPos));
    MS2map(and(PolPos > MS2_end(i), PolPos <L_MS2)) = MS2_fluorval; %Pol II's that are past the MS2 sequence.
    frac_ind = and(PolPos > MS2_start(i), PolPos < MS2_end(i)); %Indices for partially transcribed loops.
    MS2map(frac_ind) = (PolPos(frac_ind) - MS2_start(i)).*MS2_fluorval./(MS2_end(i) - MS2_start(i));

    MS2 = MS2 + sum(MS2map,2)';
    
    %Termination rate
    MS2term = zeros(size(PolPos)); %Fluorescence map containing Pol II's that are past end of gene.
    
    %For each Pol II that is part the end of gene, remove its fluorescence
    %at the termination rate. That is: F_remove = r_termination * t_past,
    %where t_past is the time it has spent past the end of gene, calculated
    %by dividing the length past the end of gene by the elongation rate.
    IndPastL = find(PolPos > L_MS2); %Indices of Pol II past end of gene.
    MS2term(IndPastL) = MS2_fluorval - r_termination * (PolPos(IndPastL) - L_MS2)./v_elong;
    MS2term(MS2term < 0) = 0; %Replace negative fluorescences with zero.
    
    %Add this fluorescence to the overall MS2 signal
    MS2 = MS2 + sum(MS2term,2)';
    
    %Basal activity
    MS2(MS2<MS2_basal) = MS2_basal;

    %Calculate PP7 signal by making a fluorescence map
    PP7_fluorval = PP7_loopn(i)/24;
    PP7map = zeros(size(PolPos));
    PP7map(and(PolPos > PP7_end(i), PolPos <L_PP7)) = PP7_fluorval; %Pol II's that are past the MS2 sequence.
    frac_ind = and(PolPos > PP7_start(i), PolPos < PP7_end(i)); %Indices for partially transcribed loops.
    PP7map(frac_ind) = (PolPos(frac_ind) - PP7_start(i)).*PP7_fluorval./(PP7_end(i) - PP7_start(i));

    PP7 = PP7 + sum(PP7map,2)';
    
    %Termination rate
    PP7term = zeros(size(PolPos)); %Fluorescence map containing Pol II's that are past end of gene.
    
    %For each Pol II that is part the end of gene, remove its fluorescence
    %at the termination rate. That is: F_remove = r_termination * t_past,
    %where t_past is the time it has spent past the end of gene, calculated
    %by dividing the length past the end of gene by the elongation rate.
    IndPastL = find(PolPos > L_PP7); %Indices of Pol II past end of gene.
    PP7term(IndPastL) = PP7_fluorval - r_termination * (PolPos(IndPastL) - L_PP7)./v_elong;
    PP7term(PP7term < 0) = 0; %Replace negative fluorescences with zero.
    
    %Add this fluorescence to the overall MS2 signal
    PP7 = PP7 + sum(PP7term,2)';
    
    %Basal activity
    PP7(PP7<PP7_basal) = PP7_basal;

end
    
end