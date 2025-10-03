function S=wheatstone_geometry(Sin)
% usage: S=wheatstone_geometry(S)
% this function generates a geometry of a vascular network, consisting of
% connectivity and length of vascular segments. 
%
% S is a model definition, this geometry function adds the geometry to the
% model. 
% This particular geometry file serves as a example. 
% 
%   IE: an array of structure, each structure is an internal elements, i.e. a segment 
%   within the network that is not directly connected to source or sink
%       nodes:  the connected nodes (horizontal vector)
%       l:      the length (m)
%       r:      the radius (m)
%       G:      the conductance (S.I.)
%       additional fields can be added in this m file (shear stress, XYZ coordinates, wall
%       thickness), depending on purpose
%       hemodynamics are not calculated here
%   SE: an array of structure, each structure us a source element, which is connected to a source (or sink) at
%   one end and internal elements on the other end
%       same fields as above
%       sourceP: e.g. 13300  the source (or sink) pressure (S.I.)
%	input parameter S: any relevant parameters, in this case fluidviscosity

% NB when using in adaptation, the r would typically be the initial value,
% or be irrelevant because the initial state is defined elsewhere (in the 5
% state, r is not a state but strain and resting radius are). Also, G would
% be calculated later anyway. Finally the input pressures could be defined
% in the default simulation parameters rather than here. So, basically,
% this file could only define the lengths, connectivity and positions of
% start and end points when used in simulations of adaptation. 

% keep the current fields in the model
S=Sin;

%% Define the version of this geometry
% note that this file contains choices and parameters, if these change and if these changes do not appear
% in S, a new version number is needed.
% new version number is need.
S.geom.version='wheatstone version aug 19 2021';

%% ==== THE CASE ===========
% wheatstone bridge case, see
% https://en.wikipedia.org/wiki/Wheatstone_bridge with the current meter
% replaced by an extra collateral-like segment
% source  is connected to internal node 1 by SE(1) 
% sink  is connected to internal node 4 by SE(2) 
% IN 1 connects to IN 2 and IN 3 by IE 2 and IE 3
% IN 2 connects horizontally to IN 3 by IE(5)
% IN 2 and IN 3 connect to IN 4 by IE(3) and IE(4)
% IN 4 connects to source  (rather a sink) of 5 mmHg via SE(2)

%% define the connectivity of the internal elements
IE(1).nodes=[1 2]; % this element connects nodes 1 and 2;  
IE(2).nodes=[1 3];
IE(3).nodes=[2 4];
IE(4).nodes=[3 4];
IE(5).nodes=[2 3]; % the collateral one

%% define how the sources and sinks are connected
% table for elements that connect to sources/sinks and a single node
SE(1).node=1;
SE(2).node=4; % note this is the venous outflow sink

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
l=1e-3;
d=sqrt(2); % keep lengths 1 in this example
dx=l*cos(pi/6); dy=l*sin(pi/6); % 60 degrees angles

IN(1).pos=[1e-3 1e-3]; % (m), x and y value
IN(2).pos=IN(1).pos + [dx dy];
IN(3).pos=IN(1).pos + [dx -dy];
IN(4).pos=IN(1).pos + [2*dx 0];

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
IE(2).r=0.9*S.r0; %  (m) <== smaller radius

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

SE(find(S.sources)).sourceP=S.sourceP;		% (N/m2)
SE(find(~S.sources)).sourceP=S.sinkP;		% (N/m2)
%% collect the elemens and nodes in S
S.IE=IE;
S.IN=IN;
S.SE=SE;

