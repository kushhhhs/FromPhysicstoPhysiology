% thrombus / thrombolysis script setup

% these are three parts:
% 1) stuff like you need to add to your project in order to be able to simulate clot
% lysis

% 2) just a quick and dirty thrombus injection loop. 

% 3) a thrombus digestion part. Every time you run this section (see Run
% Section knob in menu) it counts down. 

%% part one, assuming the 2x1 honeycomb 

% quickly tweaking the geometry, because now all vessels are 60 microns...
[S.IE.r]=vout(1E-6*[90 90 80 70 60 50 40 30 70 60 50]);

% changing the radii requires recalculating the hemodynamics
[S.IE(:).G]=vout(conductance([S.IE(:).r],[S.IE(:).l],S.fluidviscosity)); % conductance G recalculated for the three elements
[S.IN,S.IE, S.SE]=solvehemodyn(S.IN,S.IE,S.SE);

% Let's count down from 10
S.ThrombusCounterStart=10; % !! this is a model parameter, needs to go into the modelpar definitions!

% this is the counter for how long the thrombus will last. For nodes
% without a thrombus, we set it to -1;
[S.IN(:).ThrombusCounter]=deal(-1);

% we simulate the thrombus in a node by setting the radii of all connected vessels
% to almost zero, but now we need to know what they were originally.

[S.IE(:).origRadius]=vout([S.IE(:).r]); % remember the original radii

%% part 2: let's get some microthrombi in

vthrombusRadius=1e-6*[41 70 30 80]; % so four thrombi

for thr=vthrombusRadius
	[S,nodeStuck]=ThrombusFlowThrombolysis (S,thr); % use this version of ThrombusFlow as it includes the ThrombusCounter
	nodeStuck % just display which node became occlused
end


%% part 3, thrombolysis part, each run of this code counts down by 1
% this piece of code counts down for all thrombi that are currently stuck. If one
% reaches ThrombusCounter=0, it is dissolved by thrombolysis, and
% we have to remove it.

[S.IN(:).ThrombusCounter]=vout([S.IN(:).ThrombusCounter]-1); % we count down for all nodes, even if there's no thrombus, being lazy
vDissolved=find([S.IN(:).ThrombusCounter]==0); % this is a vector with the indices in IN where we need to remove the thrombus since the counter reached zero

[S.IN(:).ThrombusCounter] % print out the thrombus counter for all nodes for checking, ignore the negative countdown since these nodes are not occluded

for in=vDissolved % for each element n in vDissolved
	S.IN(in).occupied=false; % no longer occupied
	S.IN(in).ThrombusCounter=-1; % defined as no thrombus
	ie=S.IN(in).ie; % indices of connected internal elements
	[S.IE(ie).r]=vout([S.IE(ie).origRadius]); % restore the original radius of these elements
end

if ~isempty(vDissolved) % if we restored radii
	% since we changed the radii, we have to recalculate the
	% conductance and again solve flow in the whole network
	[S.IE(:).G]=vout(conductance([S.IE(:).r],[S.IE(:).l],S.fluidviscosity)); % G recalculated for the three elements
	[S.IN,S.IE, S.SE]=solvehemodyn(S.IN,S.IE,S.SE);
end

% total inflow should go up if thrombi are dissolved
 disp (['inflow (m3/s): ', num2str(S.SE(1).Q)])
