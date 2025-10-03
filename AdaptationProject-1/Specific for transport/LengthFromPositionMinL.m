function IE=LengthFromPositionMinL(IE,IN,minL)
% calculates the length of the segments based on the position of the nodes,
% works for 2D and 3D
% 1 nov 2022: added third parameter minimum length in [m]
for i=1:length(IE)
	IE(i).l=norm(IN(IE(i).nodes(1)).pos-IN(IE(i).nodes(2)).pos);
end
if nargin==3
    [IE.l]=vout(max([IE.l],minL));
end