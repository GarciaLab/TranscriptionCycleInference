function [ElongationSegments,stemloops,x_stall] = library(construct)
% Construct library. Input: construct identifier string; Output: construct
% hyperparameters

%% Define construct parameters.

if strcmp(construct,'P2P-MS2v5-LacZ-PP7v4')

    % Construct segments of potentially different elongation rate
    segments = 6.626; %Length of construct in kb
    velocities = {'v'}; %Name of elongation rate parameter

    % Stem loop positions
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 4.292; %Position of start of PP7 loops
    PP7_end = 5.758; %Position of end of PP7 loops

    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 24; %Number of PP7 loops

    % Stalling sites, where premature termination may occurr
    x_stall = []; %No stalling sites, no premature termination

    % Number of consistently bound fluorophores during the whole observation time
    MS2_0 = 0;
    PP7_0 = 0;

elseif strcmp(construct,'test')

    % Construct segments of potentially different elongation rates
    segments = [0.024,1.299,4.292,5.758,6.626]; %Segments of construct in kb (values determine ends of consecutive segments, the last entry segments(end) is always the length of the construct)
    velocities = {'v1','v2','v1','v2','v1'}; %Assignment of elongation rates to segments

    % Stem loop positions
    MS2_start = 0.024; %Position of start of MS2 loops
    MS2_end = 1.299; %Position of end of MS2 loops
    PP7_start = 4.292; %Position of start of PP7 loops
    PP7_end = 5.758; %Position of end of PP7 loops

    MS2_loopn = 24; %Number of MS2 loops
    PP7_loopn = 24; %Number of PP7 loops

    % Stalling sites, where premature termination may occurr
    x_stall = 2.900; % Sites, where premature termination (including stalling) is supposed to occur

    % Number of consistently bound fluorophores during the whole observation time
    MS2_0 = 0;
    PP7_0 = 0;

    %Define your own custom construct by modifying the above values and
    %changing the appropriate values (see template below)
    %    segments = <VECTOR OF SEGMENT ENDS>; %Segments of construct in kb (values determine ends of consecutive segments, the last entry segments(end) is always the length of the construct)
    %    velocities = <CELL ARRAY OF ELONGATION RATE PARAMETER NAMES>;
    %    %Assignment of elongation rate parameters to segments (different segments can have the same elongation rate)
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
    %    % Stalling sites, where premature termination may occurr
    %    x_stall = <VECTOR OF STALLING SITES>; % Sites, where premature termination (including stalling) is supposed to occur

    %  % Number of consistently bound fluorophores during the whole observation time
    %   MS2_0 = 0;
    %   PP7_0 = 0;

else
    error("Incorrect construct identifier.");
end

%% Output
ElongationSegments = struct('segments',{segments},'velocities',{velocities});
stemloops = struct('MS2_start',{MS2_start},'MS2_end',{MS2_end},'PP7_start',{PP7_start},'PP7_end',{PP7_end},'MS2_loopn',{MS2_loopn},'PP7_loopn',{PP7_loopn},'MS2_0',{MS2_0},'PP7_0',{PP7_0});