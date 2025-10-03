function [S,nodeStuck]=ThrombusFlow (S,thrombusRadius)
% ThrombusFlow simulates the flow of a thrombus as defined in
% thrombusRadius. The function returns the node humber where the thrombus
% got stuck (empty if the thrombus could exit the network) and updates the model S
% such that all elements connected to that node get a minimal diameter (to mimic
% flow obstruction by the clot) and recalculation of hemodynamics 


nodeStuck=[]; % return empty if thrombus didn't get stuck

% start from the node connected to the source element
currentNode=S.SE(1).node; % current node is first one
stuck=0; % not yet stuck
while ~stuck
	cn=S.IN(currentNode).cn; % all connected nodes
	ie=S.IN(currentNode).ie; % and the IE that connects them to the CurrentNode
	% for e.g. the honeycomb geometry: 
	% there are 3 nodes connected to the current node, 1-2 of them
	% are downstream, downstream nodes have lower P by definition, so
	% calculate deltaP for all connected nodes
	deltaP=[S.IN(cn).P]-S.IN(currentNode).P; % negative for downstream nodes
	% BUG 2024???
	absQ=abs([S.IE(ie).Q]); % vector of all flows, absolute, deltaP indicates direction (flow direction cannot be deduced from sign, as vessel direction is arbitrary)
	% possible branches where the thrombus might go to: downstream and
	% large enough
	possible=(deltaP<0)&([S.IE(ie).r])>thrombusRadius; % e.g. [0 1 1] means only to second and third branch
	ipossible=find(possible); % e.g. [2 3] the second and third branch

	if isempty(ipossible)
		% this means there are no downstream ie that are large enough, but
		% we might have reached a node that is connected to a downstream
		% sink element. 
		if	S.IN(currentNode).nsources>0  % there is a source, let's assume it is a sink
			if max(S.SE(S.IN(currentNode).se(:)).r)>thrombusRadius % and at least one of the connected sinks is large enough
				break % break out of the whole while loop
			end
		end
		
		% if we didn't break out above, the thrombus cannot go anywhere:
		stuck=1;
		nodeStuck=currentNode; % part of our output
		S.IN(currentNode).occupied=true; % and in the tree, keep track on which node became occupied

		% make the radii of all connected branches very small to
		% simulate their closure by the thrombus
		% this is both downstream and upstream
		[S.IE(ie).r]=deal(1e-8); % [m] 

		% since we changed the radii, we have to recalculate the
		% conductance and again solve flow in the whole network
		[S.IE(ie).G]=vout(conductance([S.IE(ie).r],[S.IE(ie).l],S.fluidviscosity)); % G recalculated for the three elements
		[S.IN,S.IE, S.SE]=solvehemodyn(S.IN,S.IE,S.SE);
		
	elseif length(ipossible)==1
		% can only go one way: go to this node
		currentNode=cn(possible);
	else
		% two possibilities, assuming 3-nodes only
		Prob=absQ(ipossible(1))/sum(absQ(possible)); % thrombi go with the flow
		% generate a choice based on these probabilities
		x=rand; % random number homogenous between 0-1
		ichoice=2-(x<Prob); % 1 or 2 for first or second possibility, e.g. 20% of flow to first downstream node means 20% chance the thrombus goes there
		currentNode=cn(ipossible(ichoice)); % and go there
	end
end
 