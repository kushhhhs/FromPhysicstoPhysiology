%% === PARAMETER SWEEP: DIFFERENT ADAPTATION SPEEDS (k) ===
clear; clc; close all;

% === Select which plots to display ===
showRadius = true;   % set to false to hide radius plots
showFlow   = false;  % set to false to hide flow plots
showWSS    = false;  % set to false to hide WSS plots

% Define the k values to test (adaptation rates)
k_values = [0.0125, 0.01, 0.05, 0.03];

% Preallocate storage for comparison
final_radii = cell(length(k_values),1);

for ik = 1:length(k_values)
    S = struct();
    S.Assignment = 2;
    warning('off','MATLAB:singularMatrix')
    warning('off','MATLAB:nearlySingularMatrix')
    S.SingularMatrixWarning = 'off';

    % --- 1. Model setup ---
    S = ModelPars(S);
    S = DefineTopology(S);

    [S.SE(find(S.sources)).sourceP] = deal(S.sourceP);
    [S.SE(find(~S.sources)).sourceP] = deal(S.sinkP);
    [S.IE.r] = vout([S.IE.r0]);
    [S.IE.G] = vout(conductance([S.IE.r],[S.IE.l],S.fluidviscosity));
    [S.IN,S.IE,S.SE] = solvehemodyn(S.IN,S.IE,S.SE);
    [S.IE.WSS] = vout(calcshearstress([S.IE.Q],[S.IE.r],S.fluidviscosity));

    % --- 2. Simulation parameters ---
    S.tend = 15;                % duration of the steady adaptation phase [s]
    S.k = k_values(ik);         % adaptation rate [1/s]
    S.WSSref = 5;               % reference wall shear stress [Pa]
    options = odeset('RelTol',1e-4,'MaxStep',1);
    X0 = [S.IE.r0]';            % initial radius vector

    % --- 3. STEADY ADAPTATION RUN ---
    [tout, mX] = ode45(@(t,X) adapt(t,X,S), [0 S.tend], X0, options);

    % Store final steady-state radii
    final_radii{ik} = mX(end,:);

    % --- 4. PULSATILE SIMULATION (starting from steady state) ---
    X0_pulse = mX(end,:);
    t_current = 0; dt = 0.01;
    pulse_amp = 0.2; pulse_freq = 1;    % amplitude [fraction], frequency [Hz]
    pulse_period = 1/pulse_freq;
    tend_pulse = 200*pulse_period;      % total simulation duration
    tout_pulse = t_current:dt:tend_pulse;

    % Preallocate arrays
    nElem = length(X0_pulse);
    mX_pulse   = zeros(length(tout_pulse), nElem);
    mQ_pulse   = zeros(length(tout_pulse), nElem);
    mWSS_pulse = zeros(length(tout_pulse), nElem);

    Xit = X0_pulse;
    it = 1;

    while t_current <= tend_pulse
        % Sinusoidal inlet pressure
        pulsatingP = S.sourceP * (1 + pulse_amp * sin(2*pi*pulse_freq*t_current));
        [S.SE(find(S.sources)).sourceP] = deal(pulsatingP);
        [S.SE(find(~S.sources)).sourceP] = deal(S.sinkP);

        % Small ODE step for adaptation
        tspan = [0 dt];
        [~, X_step] = ode45(@(t,X) adapt(t,X,S), tspan, Xit);
        Xit = X_step(end,:)';
        [S.IE.r] = vout(Xit);
        [S.IE.G] = vout(conductance(Xit,[S.IE.l],S.fluidviscosity));
        [S.IN,S.IE,S.SE] = solvehemodyn(S.IN,S.IE,S.SE);

        % Store results
        mX_pulse(it,:)   = Xit';
        mQ_pulse(it,:)   = [S.IE.Q]';
        mWSS_pulse(it,:) = vout(calcshearstress(mQ_pulse(it,:), Xit, S.fluidviscosity));

        % Advance time
        t_current = t_current + dt;
        it = it + 1;
    end

    %% --- 5. Selective plotting ---
    colors = lines(size(mX_pulse,2)); % unique color per vessel

    % === Radius plot ===
    if showRadius
        figure;
        hold on;
        for el = 1:size(mX_pulse,2)
            plot(tout_pulse, mX_pulse(:,el), 'LineWidth', 1.5, 'Color', colors(el,:));
        end
        hold off;
        xlabel('Time [s]');
        ylabel('Radius [m]');
        title(['Radius under Pulsatile Pressure (k = ' num2str(S.k) ')']);
        legend(arrayfun(@(x) sprintf('Vessel %d', x), 1:size(mX_pulse,2), 'UniformOutput', false), ...
               'Location', 'bestoutside');
        grid on;
    end

    % === Flow plot ===
    if showFlow
        figure;
        hold on;
        for el = 1:size(mQ_pulse,2)
            plot(tout_pulse, mQ_pulse(:,el), 'LineWidth', 1.5, 'Color', colors(el,:));
        end
        hold off;
        xlabel('Time [s]');
        ylabel('Flow [m^3/s]');
        title(['Flow under Pulsatile Pressure (k = ' num2str(S.k) ')']);
        legend(arrayfun(@(x) sprintf('Vessel %d', x), 1:size(mQ_pulse,2), 'UniformOutput', false), ...
               'Location', 'bestoutside');
        grid on;
    end

    % === Wall Shear Stress (WSS) plot ===
    if showWSS
        figure;
        hold on;
        for el = 1:size(mWSS_pulse,2)
            plot(tout_pulse, mWSS_pulse(:,el), 'LineWidth', 1.5, 'Color', colors(el,:));
        end
        hold off;
        xlabel('Time [s]');
        ylabel('Wall Shear Stress [Pa]');
        title(['WSS under Pulsatile Pressure (k = ' num2str(S.k) ')']);
        legend(arrayfun(@(x) sprintf('Vessel %d', x), 1:size(mWSS_pulse,2), 'UniformOutput', false), ...
               'Location', 'bestoutside');
        grid on;
    end

    disp(['Finished simulation for k = ', num2str(S.k)]);
end
