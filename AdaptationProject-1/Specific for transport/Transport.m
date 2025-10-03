function S=Transport(S,transportOrder)
% Simulates transport of a contrast agent through the circulation

%% redefine the direction of all elements such that Q is positive
% only do this for cases where flow is static!
for ie=1:S.nie
	if S.IE(ie).Q<0
		in1=S.IE(ie).nodes(1);
		in2=S.IE(ie).nodes(2);
		S.IE(ie).nodes=[in2,in1];
		S.IE(ie).Q=-S.IE(ie).Q;
	end
end
	

%% Calculate fluxes through the nodes
% calculates the fluxes of contrast through all the nodes
% work on internal elements in transport order

% NB we define fluxes here as mol/s and not as mol/(m2 x s) 

% initiate all flux vectors to zero
[S.IN.Jin]=deal(zeros(S.nt,1));

% calculate flux through the first node, in [mol/s]
S.IN(S.SE(1).node).Jin=S.Cin*S.SE(1).Q; %    

% derive for all the internal elements the distal contrast flux from the
% proximal contrast flux. Add these distal fluxes to the influx of the
% downstream node. Do this in transport order: only calculate Jdist if
% Jprox is already calculated
for k=1:S.nie
    ie=transportOrder(k);	 % select the next ie
    in1=S.IE(ie).nodes(1);	 % for brief, have in1 and in2 as the both internal nodes of this ie
    in2=S.IE(ie).nodes(2);
    Q=S.IE(ie).Q;			 % also for brief, have the Q through this element

    % the outflux from the proximal node into our element depends on the influx to that node and the
    % ratio of Q in the current element and Q towards the upstream node.
	% so if e.g. 30% of the upstream flow goes to the current node, also
	% 30% of contrast at each time goes here
	% in a diverging (arterial) branch this is <1, in a converging (venous) branch
    % or in a 2-node this is unity. The field IE.mothers has the indices
    % into IE of the upstream elements
	
    m=S.IE(ie).mothers;
    if isempty(m) % no upstream mother internal elements; this ie must be a very proximal aartery connected to an incoming SE
        Qin=S.SE(1).Q; % source element flow to start with
	else
		% add up the flows of all mothers, in an artery that is typically
		% one mother, in a vein typically two mothers
        Qin=0; 
        for i=1:length(m)
            Qin=Qin+S.IE(m(i)).Q;
		end
		
	end
	% the flux of contrast into ie depends on how much of the total flow in
	% the mothers goes into ie
    Jprox=S.IN(in1).Jin*Q/Qin; % Q is flow in current segment, Qin is all the flow towards the proximal node.

    % THIS IS THE TRANSPORT MODEL, assuming constant velocity over the
    % cross-sectional area and no diffusion, only advection 
	% calculate transport characteristics for single elements
	% the model: a simple time delay, no dispersion
	l=S.IE(ie).l; % length
	S.IE(ie).vel=S.IE(ie).Q/(pi*S.IE(ie).r^2); 
	S.IE(ie).tdelay=l/S.IE(ie).vel; % this is in seconds delay time, lenght/velocity
	
	% so now we know the delay time, but we can only deal with integer time
	% steps if we simulate this by a shift in C(t). There are more
	% sophisticated ways to do this but for small time steps and a long
	% C(t) vector it doesn't matter that much
	Ndelaysteps=round(S.IE(ie).tdelay/S.dt); % this is how far we should shift the vector
    if Ndelaysteps>S.nt % the delay is too long to be covered by the size of the C(t) vector
        Ndelaysteps=S.nt;
        disp('Transport: delay too long')
        ie % just showing in which ie this occurs
    end
    
	Jdist=[zeros(Ndelaysteps,1); Jprox(1:end-Ndelaysteps)]; % this is the actual shift of the distal flux from the proximal flux in ie 
	% if the delay is too long, the vector needs to be truncated
	Jdist=Jdist(1:S.nt);
	
	S.IN(in2).Jin=S.IN(in2).Jin+Jdist; % adds this flux vector to the possible already other influx of our distal node from elsewhere
end
for in=1:S.nin
    % calculate the area under the J curve, this should be zero if all the
    % contrast has left the node
	S.IN(in).Jarea=sum(S.IN(in).Jin)*S.dt; 
end

%% calculate tissue flows and fluxes
% this calculates the total flux into and out of each area
% create the J field in the tissue areas
[S.TA.J]=deal(zeros(S.nt,1));

% this also calculates the total inflow. 
% total flow at any time is zero per definition, what goes in goes out, no
% storage of blood. What counts is all incoming flow.
% here all incoming flow per area is summated. Dependiing on the geometry
% definition, the elemens are labeled as artery/capillary/vein and only
% incoming arterial or capillary flow would contribute to oxygen delivery,
% with outflowing arterial flowto be deducted from the oxygen delivery. This 
% oxygen delivery has not been implemented yet. 


[S.TA.Q]=deal(0); % all flows in areas set to zero, arterial inflows are then added below

for in=1:S.nin
	if ~isempty (S.IN(in).TissueArea) % this is a boundary node
		ta=S.IN(in).TissueArea; % multiple boundaries possible
		
			% add the flux vector, note this adds the whole S.nt-element or
			% so vector, no need to loop over dt
			% the sign of flux might be negative or positive
			
			
			% ignore earlier determination of RunsIn because it is topology
			% based and not flow based, figure out again here
			
			% all boundary nodes are 2-nodes, if the pressure halfway the InsideEdge is 
			% lower than in the node, the boundary node is an
			% entrance node
			Pinside=S.IE(S.IN(in).InsideEdge).P;
			S.IN(in).RunsIn=(S.IN(in).P>Pinside);
				
			
			if S.IN(in).RunsIn
				S.TA(ta).J=S.TA(ta).J+abs(S.IN(in).Jin); 
				% flow is in the IE structure, take one of the two
				% connected IE, they both have the same flow (all boundary nodes are 2-node)
				ie=S.IN(in).ie(1); 
				S.TA(ta).Q=S.TA(ta).Q+abs(S.IE(ie).Q);
			else
				S.TA(ta).J=S.TA(ta).J-abs(S.IN(in).Jin); 
            end
    end
end

% determine the concentration over time for each tissue area; THIS MISSES
% AREA SIZE!!!
[S.TA.C]=deal(zeros(S.nt,1));
for ta=1:S.N_Tissue_Areas
	S.TA(ta).C=cumsum(S.TA(ta).J)*S.dt; 
end


