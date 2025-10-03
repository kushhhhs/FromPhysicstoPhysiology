function [PLIntersection]=PolyLineIntersection(P,L)
% EVB october 2022
% finds all intersection points of line and polygon
% indicates inward or outward intersection, inward being that first point
% of line is outside polygon
% Polygon n-by-2 [x1,y1;...;xn,yn];
% Line 2-by-2 [xa,ya;xb,yb];
% PLIntersection: array of structure
% x,y: position of intersection
% iedge: index of edges crossing
% inward: array of 1/0, 1 for inward, for each crossing
nedge=size(P,1); % edge is line, vertex is node
P=[P;P(1,:)]; % copy first node to extra node to easily create the last edge
icross=0; % no crossings yet
CrossingIn=~inpolygon(L(1,1),L(1,2),P(:,1),P(:,2)); % 1 if first point of line is outside polygon; standard routine, can be pretty complex if polygon is self-intersecting
PLIntersection=[];
for iedge=1:nedge
	[xi,yi] = linexline(P(iedge:iedge+1,:),L);
	if ~isempty(xi)
		icross=icross+1;
		PLIntersection(icross).x=xi;
		PLIntersection(icross).y=yi;
		PLIntersection(icross).iedge=iedge;
		PLIntersection(icross).inward=CrossingIn;
		CrossingIn=~CrossingIn; % by definition line cannot go in or out twice in a row
	end
end


%function [xi,yi] = linexline(L1x, L1y, L2x, L2y, showIntersectionPlot)
function [xi,yi]=linexline(PL1,LL2)
% first argument part of the polygon
% Data
x1 = PL1(1,1); y1 = PL1(1,2);
x2 = PL1(2,1); y2 = PL1(2,2);
x3 = LL2(1,1); y3 = LL2(1,2);
x4 = LL2(2,1); y4 = LL2(2,2);

%------------------------------------------------------------------------------------------------------------------------
% Line segments intersect parameters
u = ((x1-x3)*(y1-y2) - (y1-y3)*(x1-x2)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));

%------------------------------------------------------------------------------------------------------------------------
% Check if intersection exists, if so then store the value
if (u >= 0 && u <= 1.0) && (t >= 0 && t <= 1.0)
    xi = ((x3 + u * (x4-x3)) + (x1 + t * (x2-x1))) / 2; 
    yi = ((y3 + u * (y4-y3)) + (y1 + t * (y2-y1))) / 2;
else
    xi = []; %NaN;
    yi = []; %NaN;
end
