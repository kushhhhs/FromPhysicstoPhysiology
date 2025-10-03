function rdot=adapt(t,r,S)
% function for ODE simulation
% t is the time r is a vector with the radii of all internal elements (i.e. the state
% variable, e.g. 5 elements for wheatstone). S is the whole model.
% rdot is a vector with the time derivatives of all radii, e.g. 5 elements

% if this function is called with a second output parameter Y, that will
% provide any derived (non-state) variables. 

%% CALCULATE REQUIRED DERIVED VARIABLES ON PER SEGMENT BASE THAT ARE NEEDED FOR NETWORK CALCULATIONS
G = conductance(r,[S.IE.l]',S.fluidviscosity);

%% CALCULATE NODE PRESSURES AND ELEMENT FLOWS, WALL STRESS, SHEAR STRESS
% This requires embedding variables for network calculations in the S.IE
% structure
[S.IE.G]=vout(G); % embed the conductances

% solvehemodyn uses the conductances and source pressures, and returns for
% all nodes the pressures as field and for all elements the flow and mean
% (midway) pressure 
[S.IN,S.IE, S.SE] = solvehemodyn(S.IN,S.IE,S.SE);
Q=[S.IE.Q]'; % flow through all IE, column vector
% P=[S.IE.P]'; % midway pressure in all IE not needed for the FP2P
% assignment
WSS=abs(calcshearstress(Q,r,S.fluidviscosity)); % 2024: forgot to take the absolute value


%% CALCULATE THE TIME DERIVATIVES OF THE STATE(S) 
% adapt the radii based on the WSS; if a radius gets smaller than 0.1
% micron, stop adapting in order to avoid G=0 and a Not a Number in solving
% the flows. 
rdot=(r>1e-7).*S.k.*r.*(WSS/S.WSSref-1); % [m/s]

end 



