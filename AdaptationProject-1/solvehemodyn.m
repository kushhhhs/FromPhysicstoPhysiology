function [IN,IE,SE]=solvehemodyn(IN,IE,SE)
%version 20 october 2016 Janina Schwarz

%% Usage
% This function solves the flow and pressure distribution throughout a network
% defined by nodes and elements.
%
% * IN: array of struct for internal nodes, requires fields _cn_ (connected nodes),
% _ie_ and _se_
% * IE: array of stuct for internal elements, requires fields _r_ (radius), _G_ (conductance) and
% _nodes_
% * SE: array of stuct for source elements, requires fields _r_ (radius), _G_ (conductance), _sourceP_ and
% _node_
%
% S.I. units are used: radius (_r_) in m, pressure (_P_) in N/m2, flow
% (_Q_) in m3/s, conductance (_G_) in S.I. units and wall shear rate (WSR) in
% 1/s.
%% Background
% At every node $j$, net inflow equals net outflow, i.e.
%
% $$\sum_{n\neq j} (P_j-P_n)*G_{jn} = 0$$
% 
% where $P_n$ is the pressure at node $n$ and
% $G_{jn}$ is the conductance of the element connecting nodes $j$ and $n$.
% All intenal nodes together yield a system of linear equations. This
% linear system can be rewritten in matrix notation with split
% contributions from connected internal nodes (at unknown pressure _P_) 
% and source elements (with source nodes at known pressure _Ps_):
% _A*P+B*Ps=0_.
%% Definition of A, B & Ps
% A is a system matrix defining internal connectivity.
% Diagonals reflect currents flowing away from the node, other elements currents
% flowing towards nodes.
%
% B is a matrix defining connectivity to sources. B has the same number of
% rows as A and one column per pressure source, each element is a
% conductance to a pressure source.
% [js]: Please note that I'm using a slightly different definition of A & B
% than in _solvehemodyn.m_.
% This definition also works for multiple source elements connected to one
% node and for multiple elements connecting the same two nodes.

nIN=length(IN); nIE=length(IE); nSE=length(SE);
if ~isfield(IN(1),'nconnect')
    for i=1:nIN, IN(i).nconnect=length(IN(i).cn); IN(i).nsources=length(IN(i).se); end
end

% find parallel IEs
IEnodes=sort(reshape([IE.nodes],[2,nIE])',2);

% define A & B as sparse matrices
cA=1;cB=1;

% all INs
for j=1:nIN
    cAs=cA;cBs=cB; % first indices for this row
    cn=unique(IN(j).cn);
    cIEs=IN(j).ie;
    IEnodesthisIN=IEnodes(cIEs,:);
    % do something similar also for SEs?
    for i=1:length(cn)
        conIE=cIEs(logical(sum(ismember(IEnodesthisIN,cn(i)),2)));
        % computation of joint parallel conductance
        G=0;
        for k=1:length(conIE), G=G+IE(conIE(k)).G; end
        vA(cA)=G;    
        colA(cA)=cn(i);
        cA=cA+1;
    end
    d=-sum(vA(cAs:cA-1));
    % then define B
    for i=1:IN(j).nsources
        vB(cB)=SE(IN(j).se(i)).G;
        colB(cB)=IN(j).se(i);
        d=d-vB(cB);
        cB=cB+1;
    end
    % diagonal term of A
    vA(cA)=d;
    colA(cA)=j;
    rowA(cAs:cA)=j;rowB(cBs:cB-1)=j; % one row for every IN
    cA=cA+1;
end

A=sparse(rowA,colA,vA,nIN,nIN);
B=sparse(rowB,colB,vB,nIN,nSE);
clear v* col* row* IEnodes % save some memory 

% vector with source pressures
PS=sparse([SE.sourceP]');

% full(A)
% full(B)
% full (PS)
% BPS=(B*PS); full(BPS)


%% Solution by matrix inversion
% Find pressures at internal nodes and in elements:
 
% spparms('spumoni',1); % info about solver
P=-A\(B*PS); % actual solution by matrix inversion
P=full(P);

% distribute over the elements
for j=1:nIN, IN(j).P=P(j); end  

% flows and mean pressures in elements
for i=1:nIE
    P1=IN(IE(i).nodes(1)).P; P2=IN(IE(i).nodes(2)).P;
    IE(i).Q=(P1-P2)*IE(i).G;
    IE(i).P=0.5*(P1+P2);    
end

% for source elements, flow is positive into the network, venous flow thus
% shows as negative
for i=1:nSE
    P1=SE(i).sourceP; P2=IN(SE(i).node).P;
    SE(i).Q=(P1-P2)*SE(i).G;
    SE(i).P=0.5*(P1+P2);
end


