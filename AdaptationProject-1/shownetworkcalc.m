% 2024: maybe this helps a bit understanding the solve hemodynamics
% function, this is slightly outdated but if you draw the network and write
% out the Kirchhoff flow equations for each node, you should see how we
% define matrices A and B and how we solve them
% shownetworkcalc is an old MATLAB script file pointing out how to solve large
% networks of resistances



% ==== THE CASE ===========
% wheatstone bridge case, source  is 100 mmHg, connected to node 1 with 
% pressure Ps1 by conductance Gs1, 
% N1 connects to N2 and N3 by G12 and G13
% N2 connects horizontally to N3
% N2 and N3 connect to N4
% N4 connects to source  (rather a sink) of Ps4=5 mmHg via a Gs4

% some numbers for internal conductance, in (mL/min) / mmHg
G12=2;
G13=4;
G23=6;
G24=3;
G34=8;

% and a non-existing path, where Gxy=0
G14=0;

%some numbers for conductances to sources, 
Gs1=1; % connects P1 to arterial inut
Gs2=0; % nothing connected here
Gs3=0;
Gs4=2; % and P4 to venous output pressure
%===============================================
%=== Set up the matrices  =====
%===============================================
%matrices are 'normal style' here, for large networks matrices are empty at
%most of the places and you would use 'sparse matrices' in matlab

%A is a system matrix defining internal connectivity
%this is A, as can be seen by writing sum(currents)=0 for each node in
%steady state
%diagonals reflect currents flowing away from the node, other elements currents
%flowing towards nodes

% 2024 so A*P (martix multiplication) provides the Kirchhoff laws applied to the four nodes, if
% there were no sources


A=[ -(G12+G13+G14), G12,            G13,            G14;
    G12,            -(G12+G23+G24), G23,            G24;
    G13,            G23,            -(G13+G23+G34), G34;
    G14,            G24,            G34,            -(G14+G24+G34)]

%B is a matrix defining connectivity to sources, B is square with
%only diagonals filled in, each element is a conductance to a pressure
%source

B=diag( [Gs1 Gs2 Gs3 Gs4],0)


% S is the input pressures connected to each of the nodes, non-connected nodes are irrelevant and set to zero 
%the input pressures for the sources, mmHg
PS=[100;
    0;
    0;
    5]
    



% (A-B)*P= -B*PS
% and then after some messing with matrix divisions:
% use the special \ (left) matrix division

P=-(A-B)\B*PS


