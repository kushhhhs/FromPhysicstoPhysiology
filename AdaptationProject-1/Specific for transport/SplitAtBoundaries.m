function [IE,nie,IN,nin]=SplitAtBoundaries(S,IE,SE,IN,TA)
% We agreed to derive the flux of contrast into several pieces (rectangles) of tissue. 
% This is the part where I will give some help, this is what I will do:
%
% I will make a function defining the brain tissue areas, as a number of 
% rectangles next to each other that cover the brain tissue. I assume your 
% network runs from top to bottom, and I will allow you to define the vertical 
% limits S.Ytop at the arterial side and S.Ybottom at the venous side. I will 
% also allow you to define the left and right limits of the whole brain, and 
% the number of tissue areas next to each other, S.Xleft and S.Xright and 
% S.N_Tissue_Areas.
%
% I will make a function that splits up all segments in your network that 
% cross a boundary, such that we get a new node at exactly the boundary. You 
% will need that later.

% For the flow, for my added nodes, I will add IN(i).TissueAreas as a vector, 
% IN(i).Qta as the total flow to this node, also a vector. IN(3).TissueAreas=[4 5] 
% and IN(3).Qta=[-2 2] means that a flow of 2 (m3/s) moves from area 4 to 5. 
% This flow does not change over time, since we have a network with steady flow. 
% If IN(i).TissueAreas is empty, that node is not a tissue boundary node. 

% The total flow into an area would be zero: what goes in goes out. What is relevant 
% is the total arterial flow, that is perfusing the tissue areas. I will add a 
% label to the boundary nodes that is again 1 2 3 for arteries capillaries and veins. 

% this function splits all segments that cross tissue boundaries, by
% introducing a 2-node at the boundary. Specifically programmed for the
% transport assignment in FP2P 2022, also can be used in FP2P 2023


nie=length(IE); nin=length(IN); 

S.InShift=0.05; % needs to go to par file, fraction inward shift of position of new boundary node

for ta=1:S.N_Tissue_Areas
	P=TA(ta).Vertices;
	OLDnie=nie;
	for ie=1:OLDnie
		in1=IE(ie).nodes(1);
		in2=IE(ie).nodes(2);
		L=[IN(in1).pos;IN(in2).pos];
		PLIntersection=PolyLineIntersection(P,L);
		if length(PLIntersection)==2
			% properties new node
			nin=nin+1;
			xm=mean([PLIntersection(1:2).x]);
			ym=mean([PLIntersection(1:2).y]);
			IN(nin).x=xm;IN(nin).y=ym;IN(nin).pos=[IN(nin).x,IN(nin).y];
			% properties existing element
			IE(ie).nodes=[in1,nin]; % old segment now from old in1 to new nin
            % properties new element
			nie=nie+1;
			IE(nie).nodes=[nin,in2]; % new segment from nin to in2
			IE(nie).r0=IE(ie).r0; 
		elseif length(PLIntersection)>=2
			disp('IE crosses area at more than two points, not covered by code');
		end
	end
end



for ta=1:S.N_Tissue_Areas
	P=TA(ta).Vertices;
	OLDnie=nie;
    for ie=1:OLDnie % don't do stuff below for newly created edges
    
		in1=IE(ie).nodes(1);
		in2=IE(ie).nodes(2);
		L=[IN(in1).pos;IN(in2).pos];
		PLIntersection=PolyLineIntersection(P,L);
        if ~isempty(PLIntersection)
            0;
		end
		% for each intersection, make a new node
		cie=ie; % in second run loop below, address a different ie
		for iis=1:length(PLIntersection)
            nin=nin+1;
            IE(cie).nodes=[in1,nin]; % old segment now from old in1 to new nin
            nie=nie+1;
			IE(nie).nodes=[nin,in2]; % new segment from nin to in2
			IE(nie).r0=IE(cie).r0; % new segment gets the same radius
            % shift the new point a bit inward
			%IN(nin).x=PLIntersection(iis).x;
            %IN(nin).y=PLIntersection(iis).y;
			if PLIntersection(iis).inward
                IN(nin).InsideEdge=nie;% the new internal element runs inside the tissue area
                IN(nin).RunsIn=1;
				% shift point a bit inward, i.e. away from first node that
				% was outside
				IN(nin).x=IN(in1).x +(1+S.InShift)*(PLIntersection(iis).x-IN(in1).x );
				IN(nin).y=IN(in1).y +(1+S.InShift)*(PLIntersection(iis).y-IN(in1).y );
			else
                IN(nin).InsideEdge=ie;% the old ie runs inside the tissue area
                IN(nin).RunsIn=0;
				% shift point a bit inward, i.e. closer to first node that
				% was inside
				%IN(nin).x=S.InShift*IN(in1).x +(1-S.InShift)*PLIntersection(iis).x;
				%IN(nin).y=S.InShift*IN(in1).y +(1-S.InShift)*PLIntersection(iis).y;
				IN(nin).x=IN(in1).x +(1-S.InShift)*(PLIntersection(iis).x-IN(in1).x );
				IN(nin).y=IN(in1).y +(1-S.InShift)*(PLIntersection(iis).y-IN(in1).y );
			end
			IN(nin).pos=[IN(nin).x,IN(nin).y];
            IN(nin).TissueArea=ta;
			% this below should not ne necessary
			if length(PLIntersection)>1
				in1=nin;
				cie=nie;% if lines crosses polygon more than once, we need to take the new point for the second run I think
			end
		end
    end
end

% remake the nodes, use the new version of MakeNodeTable that keeps
% existing fields (like position)
[IN,nin]=MakeNodeTable(IE,SE,IN);

end % of full routine










