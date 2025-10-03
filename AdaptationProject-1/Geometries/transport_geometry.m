function S=your_transport_geometry_ed(Sin)
% usage: S=your_transport_geometry(Sin)
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
%	input parameter S: any relevant parameters, e.g. fluidviscosity

% This version is specific for the transport assignment. It also creates
% tissue areas and generates extra 2-nodes at the boundary of those areas

% keep the current fields in the model
S=Sin;

%% Define the version of this geometry
% note that this file contains choices and parameters, if these change and if these changes do not appear
% in S, a new version number is needed.
% new version number is need.
S.geom.version='your transport geometry Ed version October 30 2022';

%% ==== THE CASE ===========
% => make a network yourself, this is a simple stand in made by hand.
% area covered is x=[0.5 3.5], y=[0 2] or so


% nodes connected to each internal element
% types are arterial/capillary/venous
IE(1).nodes=[1 2 ];         IE(1).type='A';
IE(2).nodes=[2 3 ];         IE(2).type='A';
IE(3).nodes=[3 4 ];         IE(3).type='C';
IE(4).nodes=[4 5 ];         IE(4).type='V';
IE(5).nodes=[2 6 ];         IE(5).type='A';
IE(6).nodes=[6 7 ];         IE(6).type='C';
IE(7).nodes=[7 11 ];        IE(7).type='V';
IE(8).nodes=[1 8 ];         IE(8).type='A';
IE(9).nodes=[8 9 ];         IE(9).type='A';
IE(10).nodes=[9 10 ];       IE(10).type='C';
IE(11).nodes=[10 11 ];      IE(11).type='V';
%IE(12).nodes=[8 12 ];       IE(12).type='A';
IE(12).nodes=[12 8 ];       IE(12).type='A'; % DEBUG: element omkeren
IE(13).nodes=[12 13 ];      IE(13).type='C';
IE(14).nodes=[13 14 ];      IE(14).type='V';
IE(15).nodes=[5 14 ];       IE(15).type='V';
IE(16).nodes=[11 5 ];       IE(15).type='V';



%% define how the sources and sinks are connected
% table for elements that connect to sources/sinks and a single node
% ADAPTED FOR TRANSPORT: SOURCE AND SINK RECONNECTED, NOT A WHEATSTONE
% ANYMORE
SE(1).node=1; 
SE(2).node=14; %note this is the venous outflow sink

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
% this is just by hand a pile of vessels
IN(1).pos=[2 1.8  ];
IN(2).pos=[1.5 1.3  ];
IN(3).pos=[0.96 1.1  ];
IN(4).pos=[0.92 0.9  ];
IN(5).pos=[1.8  ,0.5];
IN(6).pos=[ 1.9 ,1.1];
IN(7).pos=[1.8  ,0.9];
IN(8).pos=[2.5  ,1.3];
IN(9).pos=[2.3  ,1.1];
IN(10).pos=[2.3  ,0.9];
IN(11).pos=[2.1  ,0.85];
IN(12).pos=[3.5  ,1.1];%IN(12).pos=[3.1  ,1.1];
IN(13).pos=[2.9  ,0.9];
IN(14).pos=[2  ,0.2];
% annoyingly, we use both a pos vector and separate x and y..
% should be possible to do this in one statement but IN.pos(1) doesn't work
for i=1:length(IN)
	IN(i).x=IN(i).pos(1);
	IN(i).y=IN(i).pos(2);
end

% SCale the positions to realsitic values (in m)

Shrink=100; % network is in meters, so scale down by factor Shrink

[IN.x]=vout([IN.x]/Shrink);
[IN.y]=vout([IN.y]/Shrink);
for i=1:length(IN)
	IN(i).pos=[IN(i).x,IN(i).y];
end
%% define the radii 
[IE.r0]=deal(95e-6); % for now, this network needs realistic radii

%% ADDED FOR TRANSPORT some extra paramters, should eventually go to ModelPars
% these are the parameters for defining the tissue rectangles, these
% parameters cover the whole tissue from left to right
% these are split horizontally  in S.N_Tissue_Areas
% for a vertical network the tissue areas are parallel
S.Ytop=1.20/Shrink;
S.Ybottom=0.8/Shrink;
S.Xleft=0.3/Shrink;
S.Xright=4.3/Shrink;
S.N_Tissue_Areas=7;


%% ADDED FOR TRANSPORT define the tissue areas
TA=MakeTissueAreas(S);

% for calculating flow and fluxes into and out of the tissue areas, we will
% introduce extra 2-nodes at exactly the boundary

%% DEBUG draw the network before we split it at the boundaries; this can be removed
S.IE=IE;S.nie=length(S.IE);
S.IN=IN;S.nin=length(S.IN);
S.SE=SE;S.nse=length(S.SE);
S.TA=TA;S.nta=length(S.TA);
DrawNetwork(S,121); % second argument is figure number
%% ADDED FOR TRANSPORT split segments at tissue boundaries
% we split all elements that cross a boundary into two elements, such that
% the position of each element with respect to the tissue is well-defined,
% either in or out, and transport into and out of the tissue areas is also well-defined. 
% this was pretty tricky to do, might still have bugs....

[IE,nie,IN,nin]=SplitAtBoundaries(S,IE,SE,IN,TA); % S used for above parameters

%% define the lengths of the internal elements / don't change this
% calculate the lengths of the segments from the positions of the nodes,
% assuming straight segments
% might be obvious in some cases like the wheatstone, but not in others

% avoid length==0 for two superimposed nodes
minL=1e-4; % [m] this should go to ModelPars
IE=LengthFromPositionMinL(IE,IN,minL); 

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
S.TA=TA;

