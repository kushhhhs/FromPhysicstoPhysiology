function WSS=calcWSS(Q,r,fluidviscosity)
% calculates wall shear stress from flow Q, radius r and fluidviscosity

% WSS=4*fluidviscosity.*[IE.Q]./(pi*[IE.r].^3));
WSS=4*fluidviscosity.*Q./(pi*r.^3);


