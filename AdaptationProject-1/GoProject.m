% GoProject should run your complete code, generating your full set of
% data, analyses and plots. We indicate here which parts need to be
% developed, which called functions need updating by you, and we will give
% some hints. 
%% Define a structure S that will contain all parameters and results in a structured way
S=struct(); 
%% Indicate which assignment is simulated
% this is for me to test out all three assignments, you can select your
% part and delete the other two from around line 86
S.Assignment=2; % 1-2-3 for clot/adapt/transport

%% SUPPRESS WARNINGS - do not change -
% you might get networks where vessels become very small, causing nearly
% singular matrices when solving the Kirchhoff laws, this helps suppressing
% these warnings. 
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
S.SingularMatrixWarning='off';%% DEFINE YOUR MODEL PARAMETERS HERE
% You will probably have several parameters in your code; all choices for
% parameters should be made here, and should be stored in S. 
% FP2P Make/adapt ModelPars for your needs. This is your default model.  
S=ModelPars(S); % set the default parameters, this function defines the whole model, which is stored in structure S

%% CREATE THE NETWORK TOPOLOGY - do not change -
% FP2P don't change this statements, it creates the whole network by
% running one of the files in Geometries (or your own geometry file). The
% geometry file to use is defined above in ModelPars. 
S=DefineTopology(S); % set up the whole connectivity of the network, based on parameters given in above function
% S now contains anything we need for the default model, any
% changes from this or further specifications should be explicit in the parts below, but should not
% be set implicit in the called functions below. 

% FP2P we defined the sources and sinks, but their values were not yet
% included in the network, this is done below, don't change
[S.SE(find(S.sources)).sourceP]=deal(S.sourceP);
[S.SE(find(~S.sources)).sourceP]=deal(S.sinkP);
%% GET THE CURRENT RADIUS - do not change -
% FP2P in a static network, r remains r0. In a dynamic network, it
% changes. Initialize the actual radii to be equal to the initial radii
% that you have defined in your modelPars
[S.IE.r]=vout([S.IE.r0]);
% note this combined use of [ ] and vout to avoid looping over all IE%% Illustration of building blocks - do not change -
%% FP2P the code below provides a single simulation of a network, resolving
% the local pressures and flows. Depending on your assignment, you might
% have to iterate these CALCULATE steps (for thrombi), use them in an ODE function that provides
% dX=f(X), with X the state (for adaptation), or do a single simulation and then
% simulate/determine transport of contrast. 

%% CALCULATE REQUIRED DERIVED VARIABLES ON PER SEGMENT BASE THAT ARE NEEDED FOR NETWORK CALCULATIONS
% for some variables, this can be done on a per-segment basis, like
% conductance.
% for others, this requires processing the whole network, like pressure and
% flow. 

% this calculates the conductance of all elements based on their
% radius and length, and viscosity, Poiseuille law. 
[S.IE.G]=vout(conductance([S.IE.r],[S.IE.l],S.fluidviscosity));

%% CALCULATE NODE PRESSURES AND ELEMENT FLOWS
% solvehemodyn uses the conductances and source pressures, and returns for
% all nodes the pressures as field and for all elements the flow and mean
% (midway) pressure as field 
[S.IN,S.IE, S.SE]=solvehemodyn(S.IN,S.IE,S.SE);

%% CALCULATE FURTHER REQUIRED DERIVED VARIABLES THAT RELY ON THE NETWORK
% you may need e.g. wall shear stress: 
[S.IE.WSS]=vout(calcshearstress([S.IE.Q],[S.IE.r],S.fluidviscosity));

%% DRAW the network, very basic
% This shows a very basic drawing of your network. Edit this function to
% make nice informative images (and videos...) of your simulations:
% colours for e.g. flow, WSS, contrast concentration, thrombi passing by,
% line thicknesses showing diameter, you name it.  

DrawNetwork(S,1); % second argument is figure number
%% So now we have the initial network, 
% prior to the actual simulation of thrombi, adaptation, contrast, or other
% simulations. Below, we give some further code blocks and suggestions that
% are specific for thrombi / adaptation / transport

% remove the code parts that don't belong to your project

% we introduce some extra model parameters below, but these should go to
% ModelPars
%% For Thrombi
if S.Assignment==1
% Define a new field in the S.IN structure that indicates if it is occupied
% by a thrombus; % no nodes occupied yet
for in=1:S.nin,S.IN(in).occupied=false;end

% To illustrate the code: define a single thrombus with its radius:
% (this is a parameter, move it, and other parameters, to the ModelPars
% file)
S.thrombusRadius=50e-6; % [m]

% And let it flow into the network. The actual simulation of this is preprogrammed but
% you should understand the algorithm. You are free to change the algorithm
% if you see reason for this.
[S,nodeStuck]=ThrombusFlow (S,S.thrombusRadius); 


%% For Adaptation
elseif S.Assignment==2
% here add the adaptation, using the Matlab ODE45 function. The actual
% model should be in a odefun that produces dX/dt=adapt(t,X). 
% X here is a vector of all the radii of the internal elements, i.e. your state variable, starting at X0. 

% have some extra model parameters, we define them here for now but they
% belong in the ModelPars m file, move them to there. 
S.tend=10; % [s] end time of simulation
S.k=1; % [1/s] speed of adaptation
S.WSSref=5; % [n/m2] reference wall shear stress

% set up the simulation options, these are also model parameters!
options = odeset('RelTol',1e-4,'MaxStep',1);

% the ode45 cannot handle the S structure, you need to collect the initial radii from 
% the S.IE structure and combine them in a column vector: 
X0=[S.IE.r0]'; % column vector of initial radii

% run the simulation. here, adapt generates dX/dt based on the time and X. 
% the function is there but you need to understand it and you may want to
% adapt the adapt function, if you want to run another adaptation model
[tout,mX] = ode45(@(t,X) adapt(t,X,S), [0 S.tend], X0, options);  % <== define solver and model

% tout is a column vector with all the simulated timepoints and mX a matrix with all the radii of the ie
% at all the timepoints. The last row is the last simulation, where the system might be in steady state.  

% mX only contains the state variables. Other variables are a direct
% algebraic function of these, and can be recalculated. yet the
% recalculation of pressures and flow requires again the S.IE structure,
% this is how you could do this for derived variables in the internal elements:

mWSS=zeros(size(mX)); % predefine; r has nie columns and each row is a timestep
mQ=zeros(size(mX)); % predefine
mG=zeros(size(mX)); % predefine
for it=1:length(tout)
	Xit=mX(it,:); % radii at a given timepoint
	[S.IE.r]=vout(Xit); % let's put them back into the network
	mG(it,:) = conductance(Xit,[S.IE.l],S.fluidviscosity); % get the conductances at this timepoint
	[S.IE.G]=vout(mG(it,:)); % embed these in the network 
	[S.IN,S.IE, S.SE] = solvehemodyn(S.IN,S.IE,S.SE); % solve the network
	mQ(it,:)=[S.IE.Q]'; % extract the flows from the network
	mWSS(it,:)=calcshearstress(mQ(it,:),Xit,S.fluidviscosity);
	[S.IE.WSS]=vout(mWSS(it,:)); % embed these in the network 
	% at this point you have the network showing the state and derived variables at this
	% timepoint.   
end


%% For Contrast Transport
else % S.Assignment==3
% define the arterial input function here 
% this now is a simple step function. You may want to consider more
% realistic functions

% parameters for contrast transport, move parameter definitions to the
% ModelPars function, where they belong. 
S.dt=0.01;                     % (s) time step for tracer dilution curves
S.nt=300;                      % number of time points in curve

S.Cin=zeros(S.nt,1); 
S.Cin(1:10)=1e11; % [mol/m3]; 10^11, high because Q is very low, arterial input function concentration

% figure out the order of calculating the concentration at the nodes:
% this function you can use
[S,transportOrder]=GetTransportOrder(S);

% TransportOrder now is a vector of the internal elements 

% Now simulate the transport. This function is there for you to use but you
% need to understand the algorithms and you may want to adapt these,
% depending on your research question. 
S=Transport(S,transportOrder);

end % Assigment switch














