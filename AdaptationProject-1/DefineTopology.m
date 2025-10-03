function S=DefineTopology(Sin)
% This function sets up the topology and adds it to the model S.
% In the topology, elements are defined as vessel segments between two
% nodes (internal) or between a node and a pressure source/sink (source
% elements). Internal nodes connect internal and/or source elements. 
%
% Sin is a structure with fields that define the generation of the topology
% 
% S will be a structure with
%   
%   IE: an array of structure, each element an internal elements (i.e.
%   vessel segment), with minimally these fields:
%       nodes:  the connected nodes (horizontal vector)
%       l:      the length (m)
%
%   SE: an array of structure with source elements
%       same basic fields as above, plus
%       sourceP: the source pressure
%   IN: an array of structures with internal nodes
%       nconnect        number of connected internal elements
%       nsources        number of connected source elements
%       cn              array of connected node indices 
%       ie              array of connected internal elements
%       se              array of connected source elements
%	other fields: any network-wide parameters that define the model 

% Keep the current partially defined model
S=Sin; 

%% DEFINITION OF THE CONNECTIVITY OF THE MODEL: the element tables / DO NOT CHANGE
% LOAD THE MODEL, use the default topology if no input is defined 
if ~isfield(S,'geometry')
	S.geometry='wheatstone_geometry'; disp('using wheatstone_geometry as default model'); 
end

% fill the IE SE and IN fields as defined in the model
S=eval([S.geometry '(S)']); % user defined connectivity model, can be replaced by other models

% determine the sizes for later quick reference
S.nie=length(S.IE); % number of internal elements (segments), each connecting two nodes
S.nse=length(S.SE); % number of connections to pressure sources/sinks from single node
S.nin=length(S.IN); % number of nodes in the model, each node connects 2 or more internal or source elements

 


