function [S] = honeycomb_geometry(Sin)
% usage: S=honeycomb(S)
% Honecomb generates a honeycomb-like geometry 
% Further explanation is in the wheatstone_geometry file

% keep the current fields in the model
S=Sin;



%% Define the version of this geometry
% note that this file contains choices and parameters, if these change and if these changes do not appear
% in S, a new version number is needed.
% new version number is need.
S.geom.version='honeycomb version sep 10 2021';
%% === THE CASE ==============
% Honeycomb structure consisting of hexagonals, length and width as defined
% in the parameters hexx and hexy

% if we haven't defined this yet, define the size and number of inflow and
% outflow vessels
if ~isfield(S,'ncombsx')
	S.ncombsx  = 3; 
	S.ncombsy  = 1;
end

% define derived size indicators to get the configuration right
hexx=2+2*S.ncombsx;
hexy=S.ncombsy; 


%% define Internal Elements
for i = 3:1:hexx
    IE(i-2).nodes = [i-2, i-1];   % 1~(hexx-1)
    if hexy ==1  % if only one hexagon in y direction
        IE(i-3+hexx-1).nodes = [i-3+hexx, i-3+hexx+1];
        for k = 0:1:(hexx/2)-1
            IE(hexx*hexy + hexx -3 + k).nodes = [2*k+1, hexx+(2*k)]; % vertical element
        end
    else % more than one hexagonal
        for n = 1:1:hexy % y direction, parallel element
            IE(i-3+(n*(hexx-1))).nodes=[i-3+(n*(hexx)),i-3+n*(hexx)+1]; % nth (hexx-1)~(n+1)th (hexx-1)
            if n ~= hexy %for shift element
                IE((hexx-2)*(n+1)+n).nodes = [(hexx-1)*(n+1)+n-1, (hexx-1)*(n+1)+n];
            end
        end
        for j = 0:1:hexx/2-1 % vertical element
            for k = 1:1:hexy
                if hexy ==2
                    IE(j+1+ 1/2*hexx*(k-1) +(hexx -2)*hexy + hexx-1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2];
                else
                    if k ==1
                        IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2];
                    elseif k ==2
                        IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2+1];
                    elseif mod(k,2)==0 %even
                        if k == hexy
                            IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2];
                        else
                            IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+1+j*2, hexx*k+j*2+1];
                        end
                    else %odd
                        IE((hexx-2)*2+(hexy-1)*(hexx-1)+ hexx/2*(k-1)+j+1).nodes = [hexx*(k-1)+j*2, hexx*k+j*2];
                    end
                end
            end
        end
    end
end

%% define node connection table part one - the internal elements
% we normally would use the statement [IN,nin]=MakeNodeTable(IE,SE);
% but we now first need the nodes before we hook on the source elements
% so this is the first part of MakeNodeTable
nie=length(IE);

% total number of internal nodes
nin=max([IE.nodes]);  

%no connections at start
for i=1:nin 
    IN(i).nconnect=0;
end

% fill the internal connections based on the definition of IE,  
for ie=1:nie
    node(1)=IE(ie).nodes(1); 
    node(2)=IE(ie).nodes(2);
    for j=1:2
        k=node(j);
        IN(k).nconnect=IN(k).nconnect+1;
        nc=IN(k).nconnect;
        IN(k).cn(nc)=node(3-j); % 3-j means 1=>2 and 2=>1
        IN(k).ie(nc)=ie; % points back to element table
    end
end


%% define the position of the nodes
% in this model, the hexagons are regular, all IE have the same length,
% taken unity here and scaled at end
for i=0:1:((hexx)/2)-1              %i=0,1,2,3
    IN(1+2*i).x=2*i;                  %1,3,5,7
    IN(1+2*i).y=1;                    %1,3,5,7
    for j=0:1:((hexx)/2)-2          %j=0,1,2
        IN(2+2*j).x=1+2*j;            %2,4,6
        IN(2+2*j).y=0;                %2,4,6
        for k=0:1:hexy-1            %k=0,1
            if k==0
                IN(hexx+2*i).x=2*i;                         %8,10,12,14
                IN(hexx+2*i).y=1+sqrt(2);
                IN(hexx+1+2*j).x=1+2*j;                     %9,11,13
                IN(hexx+1+2*j).y=2+sqrt(2);
            elseif k==1
                IN(k*hexx+2*i).x=2*i;                       %8,10,12,14
                IN(k*hexx+2*i).y=k+k*sqrt(2);
                IN(k*hexx+2*i+1).x=2*i+1;                   %9,11,13,15
                IN(k*hexx+2*i+1).y=k+1+k*sqrt(2);
                IN(hexy*hexx+2*i).x=2*i+1;                %16,18,20,22
                IN(hexy*hexx+2*i).y=hexy+2*sqrt(2);
                IN(hexy*hexx+2*j+1).x=2*j+2;              %17,19,21
                IN(hexy*hexx+2*j+1).y=hexy+1+2*sqrt(2);
            else
                if mod(k,2) == 0   %number is even
                    IN(k*hexx+2*i).x=2*i;              %16,18,20,22
                    IN(k*hexx+2*i).y=k+1+k*sqrt(2);
                    IN(k*hexx+2*i+1).x=1+2*i;          %17,19,21,23
                    IN(k*hexx+2*i+1).y=k+k*sqrt(2);
                else %number is odd, so kk=3
                    IN(k*hexx+2*i).x=2*i;              %8,10,12,14
                    IN(k*hexx+2*i).y=k+k*sqrt(2);
                    IN(k*hexx+2*i+1).x=1+2*i;          %9,11,13,15
                    IN(k*hexx+2*i+1).y=k+1+k*sqrt(2);
                end
                if mod(hexy,2) == 0
                    IN(hexy*hexx+2*i).x=1+2*i;       %24,26,28,30
                    IN(hexy*hexx+2*i).y=hexy+hexy*sqrt(2);
                    IN(hexy*hexx+2*j+1).x=1+2*j+1;   %25,27,29
                    IN(hexy*hexx+2*j+1).y=hexy+1+hexy*sqrt(2);
                else
                    IN(hexy*hexx+2*i).x=2*i;         %24,26,28,30
                    IN(hexy*hexx+2*i).y=hexy+hexy*sqrt(2);
                    IN(hexy*hexx+2*j+1).x=2*j+1;     %25,27,29
                    IN(hexy*hexx+2*j+1).y=hexy+1+hexy*sqrt(2);
                end
            end
        end
    end
end
%scaling
ElementLength =1e-3; % length of each element <=== HARD_CODED MODELING CHOICE!
ratio=ElementLength/sqrt(2);

for i=1:nin
    IN(i).x=IN(i).x*ratio;
    IN(i).y=IN(i).y*ratio;
end

% current version of geometries uses IN(.).pos rather than .x and .y, such
% that we can do 3D later. 
for i=1:nin
	IN(i).pos=[IN(i).x IN(i).y]; % maybe rather have this as column vector to be able to use [IN.pos]
end
% But keep .x and .y too, or solve some indexing problems 
% rmfield(IN,{'x','y'}); % remove the x and y fields

%% Define the (initial) radii
[IE.r0]=deal(S.r0);
%% Define the source elements
% we will connect the high P source(s) to the bottom nodes, and the low P
% sink(s) to the top nodes in this model
% S.ConnectedBottomNodes is a vector of nodes to be connected to a high
% pressure source, counting simply from 1..ncombsx
% same for S.ConnectedTopNodes, also 1..ncombsx

% make a simple diagonal connection if not defined
if ~isfield(S,'ConnectedBottomNodes')
	S.ConnectedBottomNodes=1;
	S.ConnectedTopNodes=S.ncombsx;
end

% we need to find the indices in IN that are the bottom nodes and top nodes, 
% the bottom node indices into IN are simply 2*(1:S.ncombsx)
% the top node indices are nin+1-2*(S.ncombsx:-1:1); 

% (this was used to sort this out)
% S.ncombsx=8;S.ncombsy=5; S=honeycomb_geometry(S); yy=[S.IN.y];maxy=max(yy); find(yy>=(maxy-1e-31))

nbottom=length(S.ConnectedBottomNodes);
ntop=length(S.ConnectedTopNodes);

for se=1:nbottom
	SE(se).node=2*S.ConnectedBottomNodes(se);
end

% indices into IN of all top nodes:
v=nin+1-2*(S.ncombsx:-1:1);
for se=1:ntop % and find the indices in v we wabt to connect
	SE(se+nbottom).node=v(S.ConnectedTopNodes(se));
end

%% define the lengths and other properties of source-connecting elements
% we need to define length and radius and calculate conductance
% or alternatively just enter conductance
[SE.l]=deal(ElementLength); % L was defined above for determining positions, we will use it for the source element lengths too 
[SE.r]=deal(10*S.r0); % all source elements have the same radius as the starting radius of the IE 
%disp('honeycomp regel 207 aangepast 2024')
[SE.G]=vout(conductance([SE.r],[SE.l],S.fluidviscosity));

%% define node connection table part two - the source elements
nse=nbottom+ntop; 
% inflow and outflow connections are both 'source' but above are coupled to
% resp high and low pressure. For setting up the structure below this is
% not relevant
%no connections at start
for i=1:nin 
    IN(i).nsources=0;
end

% fill the connections to sources
for se=1:nse
    k=SE(se).node;
    IN(k).nsources=IN(k).nsources+1;
    IN(k).se(IN(k).nsources)=se;
end
%% define the lengths of the internal elements / don't change this
% calculate the lengths of the segments from the positions of the nodes,
% assuming straight segments
% might be obvious in some cases like the wheatstone, but not in others

IE=LengthFromPosition(IE,IN); 
%% define the external pressures
% note that these values might be changing in time, if they are chosen as
% input variables. In that case it should not be necessary to initialize
% them here. 

% the SE connected to bottom nodes get the high pressure, those on top the
% low pressures
[SE(1:nbottom).sourceP]=deal(S.sourceP);	% (N/m2)
[SE(nbottom+1:nbottom+ntop).sourceP]=deal(S.sinkP);		% (N/m2)


%% assign to the output structure
S.IE=IE; S.SE=SE;S.IN=IN;
S.nin=nin; S.nie=length(IE); S.nse=nse;


end