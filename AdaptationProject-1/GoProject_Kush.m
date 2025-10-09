S = struct();
S.Assignment = 2;
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')
S.SingularMatrixWarning='off';
% Set the default parameters of the model below

S = ModelPars(S); % Get default parameters 
S=DefineTopology(S);

[S.SE(find(S.sources)).sourceP]=deal(S.sourceP);
[S.SE(find(~S.sources)).sourceP]=deal(S.sinkP); 
[S.IE.r]=vout([S.IE.r0]); % Curretn radius fo the model 
[S.IE.G]=vout(conductance([S.IE.r],[S.IE.l],S.fluidviscosity))
[S.IN,S.IE, S.SE]=solvehemodyn(S.IN,S.IE,S.SE);% pressure at all points 
[S.IE.WSS]=vout(calcshearstress([S.IE.Q],[S.IE.r],S.fluidviscosity));% wss 
S.tend=15; % [s] end time of simulation
S.k=0.25; % [1/s] speed of adaptation
S.WSSref=5; % [n/m2] reference wall shear stress
% set up the simulation options, these are also model parameters!
options = odeset('RelTol',1e-4,'MaxStep',1);
[tout,mX] = ode45(@(t,X) adapt(t,X,S), [0 S.tend], X0, options); % define the model and the solver 
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
% Radii evolution plot
figure;
plot(tout, mX, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Radius [m]');
title('Radius Evolution');
legendStrings = arrayfun(@(x) sprintf('Element %d', x), 1:size(mX,2), 'UniformOutput', false);
legend(legendStrings);

% Wall Shear Stress evolution plot
figure;
plot(tout, mWSS, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Wall Shear Stress [N/m^2]');
title('WSS Evolution');
legend(legendStrings);

% Flow evolution plot
figure;
plot(tout, mQ, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Flow [m^3/s]');
title('Flow Evolution');
legend(legendStrings);

X0_pulse = mX(end, :);
t0_pulse = tout(end);
pulse_amp = 0.2;          % Amplitude fraction (e.g., 0.2 means ±20%)
pulse_freq = 1;           % Frequency [Hz]
pulse_period = 1 / pulse_freq;
tend_pulse = 5 * pulse_period; % simulate several cycles
dt = 0.01;                 % time step for loop
t_current = 0;
it = 1;
tout_pulse = t_current:dt:tend_pulse;

mX_pulse   = zeros(length(tout_pulse), length(X0_pulse));
mWSS_pulse = zeros(size(mX_pulse));
mQ_pulse   = zeros(size(mX_pulse));

% For full conductance matrices
nElements = length(X0_pulse);
mG_pulse   = zeros(nElements, nElements, length(tout_pulse));

Xit = X0_pulse;  % start from steady state
it = 1;
t_current = 0;

while t_current <= tend_pulse
    % Step 1: Pulsating pressure
    pulsatingP = S.sourceP * (1 + pulse_amp * sin(2*pi*pulse_freq*t_current));
    [S.SE(find(S.sources)).sourceP] = deal(pulsatingP);
    [S.SE(find(~S.sources)).sourceP] = deal(S.sinkP);

    % ODE step for radii adaptation
    tspan = [0 dt];
    [t_step, X_step] = ode45(@(t,X) adapt(t,X,S), tspan, Xit);
    Xit = X_step(end,:)';  % update radii

    % Update network
    [S.IE.r] = vout(Xit);
    Gmat = conductance(Xit,[S.IE.l],S.fluidviscosity);
    mG_pulse(:,:,it) = Gmat; % Store full conductance matrix
    [S.IE.G] = vout(diag(Gmat)); % keep diagonal for network solver
    [S.IN,S.IE,S.SE] = solvehemodyn(S.IN,S.IE,S.SE);

    %  Store results
    mX_pulse(it,:)   = Xit';
    mQ_pulse(it,:)   = [S.IE.Q]';
    mWSS_pulse(it,:) = diag(calcshearstress(mQ_pulse(it,:),Xit,S.fluidviscosity));


    t_current = t_current + dt;
    it = it + 1;
end

figure;
plot(tout_pulse, mX_pulse, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Radius [m]');
title('Radius under Pulsatile Pressure');
legendStrings = arrayfun(@(x) sprintf('Element %d', x), 1:size(mX_pulse,2), 'UniformOutput', false);
legend(legendStrings);

figure;
plot(tout_pulse, mWSS_pulse, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Wall Shear Stress [N/m^2]');
title('WSS under Pulsatile Pressure');
legend(legendStrings);

figure;
plot(tout_pulse, mQ_pulse, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Flow [m^3/s]');
title('Flow under Pulsatile Pressure');
legend(legendStrings);
fprintf('\nRadii at final time step (t = %.3f s):\n', tout_pulse(end));
for elem = 1:size(mX_pulse,2)
    fprintf('Element %d: %.6f m\n', elem, mX_pulse(end,elem));
end