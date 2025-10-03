function [S]=randomnetwork(Sin,initial)
% THIS MIGHT STILL NEED SOME DEBUGGING!

% keep existing fields
S=Sin;

%% ==== THE CASE ===========
%Generates a random network with S.nie internal elements, S.nin nodes and a single source and
%sink. Second argument is seed for random number. Network may have
%unconnected nodes, many connections etc, also connections from a node to
%itself, which don't enter the matrix equations for conductance
if nargin>2,
    rand('twister',initial); % define the random sequence
end

%Generate integers uniform on the set 1:n.
%           r = ceil(n.*rand(100,1));

for ie=1:S.nie
    IE(ie).nodes=[ceil(S.nin*rand), ceil(S.nin*rand)];
end

% table for elements that connect to sources/sinks and a single node
SE(1).node=IE(1).nodes(1);
SE(2).node=IE(S.nie).nodes(2); % could have been any but this makes sure there is a connection to inlet and outlet

% define the lengths and radii of internal elements
% and calculate conductance

% have all connections the same for now
etta=1e-3; % viscosity in SI
for i=1:S.nie
    %IE(i).l=1e-3; % (m) 
    %IE(i).r=1e-4; %  (m)
    %IE(i).G=pi*IE(i).r^4/(8*etta*IE(i).l); % conductance
    IE(i).G=10;
end

% define the lengths and other properties of source-connecting elements,
%SE(1).l=1e-3;       % m length
%SE(1).r=1e-4;      % m internal radius
SE(1).G=10; %pi*SE(1).r^4/(8*etta*SE(1).l);

%SE(2).l=1e-3;       % m length
%SE(2).r=1e-4;      % m internal radius
SE(2).G=10; %pi*SE(2).r^4/(8*etta*SE(2).l);

% define the external pressures
% define source and sink pressures in N/m2 
SE(1).sourceP=100*133; 
SE(2).sourceP=0*133; 
S.SE=SE;
S.IE=IE;
%% generate the node table / don't change this
% we now have a table of elements; we also need to have a table of internal nodes 
% this function makes a table of nodes, each node has information on the
% connected segments and will have pressure data
% never change this!
[IN,nin]=MakeNodeTable(IE,SE);
% note that the connectivity is assumed constant. If adaptation would lead
% to loss of segments or new segments, this part becomes part of the
% iterations, but since the total number of state variables changes, this
% is a major coding issue
%% define the position of the nodes
% the positions are used to calculate the lengths of the segments and are needed for
% visualization of the simulation results.
for i=1:nin
	IN(i).pos=rand(1,2); % arbitrary
end
	
	% annoyingly, we use both a pos vector and separate x and y..
% should be possible to do this in one statement but IN.pos(1) doesn't work
for i=1:length(IN)
	IN(i).x=IN(i).pos(1);
	IN(i).y=IN(i).pos(2);
end


%% define the lengths of the internal elements / don't change this
% calculate the lengths of the segments from the positions of the nodes,
% assuming straight segments
% might be obvious in some cases like the wheatstone, but not in others

IE=LengthFromPosition(IE,IN); 
%% define the (initial) radius 
[IE.r0]=deal(S.r0);

%% define the lengths and other properties of source-connecting elements
% for the simulation, only the conductance is relevant. For plotting, we
% still need to add positions. 
SE(1).l=1e-3;       % m length
SE(1).r=1e-4;      % m internal radius
SE(1).G=pi*SE(1).r^4/(8*S.fluidviscosity*SE(1).l);

SE(2).l=1e-3;       % m length
SE(2).r=1e-4;      % m internal radius
SE(2).G=pi*SE(2).r^4/(8*S.fluidviscosity*SE(2).l);

%% define the external pressures
% define source and sink pressures in N/m2 
% % SE(1).sourceP=100*133;		% (N/m2)
% % SE(2).sourceP=5*133;		% (N/m2)

% note that these values might be changing in time, if they are chosen as
% input variables. In that case it should not be necessary to initialize
% them here. 

% above we defines source at 1 and sink at 2
SE(1).sourceP=S.sourceP;		% (N/m2)
SE(2).sourceP=S.sinkP;		% (N/m2)
% SE(find(S.sources)).sourceP=S.sourceP;		% (N/m2)
% SE(find(~S.sources)).sourceP=S.sinkP;		% (N/m2)
%% collect the elemens and nodes in S
S.IE=IE;
S.IN=IN;
S.SE=SE;



