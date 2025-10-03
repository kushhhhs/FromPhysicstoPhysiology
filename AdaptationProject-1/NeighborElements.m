function [vl,vr]=NeighborElements(S,i)
% returns two vectors with connected elements to the input element S.IE(i),
% i is scalar
% source elements are not included
nodes=S.IE(i).nodes;
vl=S.IN(nodes(1)).ie;vl=vl(vl~=i);
vr=S.IN(nodes(2)).ie;vr=vr(vr~=i);

