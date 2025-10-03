function S=ModelPars(S)
% usage: S=ModelPars(S) or S=ModelPars
% This defines the default model. Without input, the mode is built from
% scratch. With input, any current fields in S are preserved.
% 
% use S as input to redefine the default without losing and having to
% rebuild the network topology 
%
% All default parameters, modeling choices, specific values should go into
% this function or in the definition of the topology. Changing the model from 
% default is only done in the GoProject main script. 

%================================================
% =========== model definition ==================
%================================================       

% redefine the model version whenever anything relevant in the code
% changes, in the main script or the function it calls
% we need to be 100% sure which code is executed. Document the changes and keep
% the old version. 

if ~nargin, S=struct(); end % use S as input to redefine the default without losing the network

S.model='FP2P model XXXX version YYYY';
S.timestamp=datestr(now);

%% model parameters
S.fluidviscosity = 0.0035;			% (Ns/m)  fluid viscosity of blood
S.sourceP= 100 * 133.322368;		% (N/m2) Input pressure to all sources, the number 133.322368 is the conversion from mmHg to N/m2
S.sinkP= 0 * 133.322368;			% sink pressures, see below

% FP2P ADD ANY OTHER PARAMETERS HERE


% FP2P choose your geometry or define a new one
% select the geometry 
% S.geometry='wheatstone_geometry'; % the geometry used for this simulation, this string will be evaluated later. See notes in the wheatstone_geometry example
if S.Assignment==3
	S.geometry='transport_geometry';
else
	S.geometry='honeycomb_geometry';
end
%S.geometry='randomnetwork_geometry';


%% Below, extra model pars can be defined that depend on the geometry


if strcmp(S.geometry,'wheatstone_geometry')
%% geometry definition for wheatstone
% anything that is specific for a geometry goes here
S.geometry='wheatstone_geometry';  % the geometry used for this simulation, this string will be evaluated later. See notes in the wheatstone_geometry example
%S.geometryversion='wheatstone version aug 19 2021'; % make sure this matches the actual geometry file version ???? why do we need this?
% have any parameters that are specific for the geometry here. 

% The geometry m file should define connectivity and fixed geometric parameters like 
% length and coordinates of elements. 

% the geometry has source elements that act indeed as sources (inflow) or rather sinks (outflow). The network solver does not discriminate between these
% but the pressures need to be defined here, and these are high for inflow
% and low for outflow. sourceP and sinkP were defined above.

% (maybe better to define the sources and sinks in the geometry)
S.r0=80e-6;							% default radius, one of the segments will be made smaller to tip the bridge out of equilibrium
S.sources=[1 0];					% indicates which of the source elements are sources or sinks, in order to set the pressures below 
									% this can only be done if the number
									% of source and sink elements is known upfront


elseif strcmp(S.geometry,'honeycomb_geometry')
%% geometry definition for honeycomb


% this geometry requires further specification here
S.ncombsx  = 2; % number of columns of hexagons, in x direction
S.ncombsy  = 1; % number of rows of hexagons

% we will connect the high P source(s) to the bottom nodes, and the low P
% sink(s) to the top nodes in this model
% S.ConnectedBottomNodes is a vector of nodes to be connected to a high
% pressure source, counting simply from 1..ncombsx
% same for S.ConnectedTopNodes, also 1..ncombsx

S.ConnectedBottomNodes=[1 ];  % so a single high P source at the bottom
S.ConnectedTopNodes=[1 2 ]; % ... and a single exit paths

% define an intial radius
S.r0= 60e-6; % FP2Pcould be radius of all vessels or max radius in network, depending on purpose 

% indicate which of the source elements are sources or sinks, in order to set the pressures below
S.sources=[ ones(1,length(S.ConnectedBottomNodes)) zeros(1,length(S.ConnectedTopNodes))]; 
elseif strcmp(S.geometry,'transport_geometry')
	S.sources=[1 0]; % one arterial inflow, one venous outflow
elseif strcmp(S.geometry,'randomnetwork_geometry')
%% geometry definition for random network
S.nin=15;
S.nie=25;
S.nse=2;
S.r0=80e-6;		
S.sources=[1 0];
%S.sinks=2;
end % geometry definition

%% For the transport function you may want to define tissue areas (in 2D)







