function G=conductance(ri,l,etta)

% calculates conductance, based on the radius and length

% parameters in structure S that are used:
% S.fluidviscosity

G=pi*ri.^4./(8*etta*l);

end