%% === PARAMETER SWEEP: DIFFERENT ADAPTATION SPEEDS (k) ===
clear; clc; close all;

% === Select which plots to display ===
showRadius = true;   % set to false to hide radius plots
showFlow   = false;  % set to false to hide flow plots
showWSS    = false;  % set to false to hide WSS plots

% === Define adaptation rate values ===
k_values = [0.001, 0.0025, 0.005, 0.01];

% === Preallocate storage ===
final_radii = cell(length(k_values),1);
WSS_data = struct();   % store time + mean WSS for each k

% === Setup figures before loop ===
if showRadius
    fig_radius = figure;
    hold on;
    colors_radius = lines(length(k_values));
end

if showFlow
    fig_flow = figure;
    hold on;
    colors_flow = lines(length(k_values));
end

if showWSS
    fig_wss = figure;
    hold on;
    colors_wss = lines(length(k_values));
end

% === Combined WSS figure setup ===
figure;
hold on;
colors = lines(length(k_values));

for ik = 1:length(k_values)
    %% --- 1. Model setup ---

    S = struct();
    S.Assignment = 2;
    warning('off','MATLAB:singularMatrix')
    warning('off','MATLAB:nearlySingularMatrix')
    S.SingularMatrixWarning = 'off';
    
    S = ModelPars(S);
    S = DefineTopology(S);
    [S.SE(find(S.sources)).sourceP] = deal(S.sourceP);
    [S.SE(find(~S.sources)).sourceP] = deal(S.sinkP);
    [S.IE.r] = vout([S.IE.r0]);
    [S.IE.G] = vout(conductance([S.IE.r],[S.IE.l],S.fluidviscosity));
    [S.IN,S.IE,S.SE] = solvehemodyn(S.IN,S.IE,S.SE);
    [S.IE.WSS] = vout(calcshearstress([S.IE.Q],[S.IE.r],S.fluidviscosity));
    WSS_before_pulse = [S.IE.WSS];
    
    % Display
    disp('💡 WSS just before pulsatile flow:');
    disp(WSS_before_pulse);
    %% --- 2. Simulation parameters ---
    S.tend = 400;                % steady adaptation duration [s]
    S.k = k_values(ik);          % adaptation rate [1/s]
    S.WSSref = 5;                % reference WSS [Pa]
    options = odeset('RelTol',1e-4,'MaxStep',1);
    X0 = [S.IE.r0]';             % initial radius vector

    %% --- 3. STEADY ADAPTATION RUN ---
    [tout, mX] = ode45(@(t,X) adapt(t,X,S), [0 S.tend], X0, options);
    final_radii{ik} = mX(end,:);
    
    %% --- 4. PULSATILE SIMULATION ---
    X0_pulse = mX(end,:);
    t_current = 0; dt = 0.01;
    pulse_amp = 0.2; pulse_freq = 1;
    pulse_period = 1/pulse_freq;
    tend_pulse = 400*pulse_period;
    tout_pulse = t_current:dt:tend_pulse;

    nElem = length(X0_pulse);
    mX_pulse   = zeros(length(tout_pulse), nElem);
    mQ_pulse   = zeros(length(tout_pulse), nElem);
    mWSS_pulse = zeros(length(tout_pulse), nElem);
    Xit = X0_pulse; it = 1;

    while t_current <= tend_pulse
        % Sinusoidal inlet pressure
        pulsatingP = S.sourceP * (1 + pulse_amp * sin(2*pi*pulse_freq*t_current));
        [S.SE(find(S.sources)).sourceP] = deal(pulsatingP);
        [S.SE(find(~S.sources)).sourceP] = deal(S.sinkP);

        % ODE step for adaptation
        tspan = [0 dt];
        [~, X_step] = ode45(@(t,X) adapt(t,X,S), tspan, Xit);
        Xit = X_step(end,:)';
        [S.IE.r] = vout(Xit);
        [S.IE.G] = vout(conductance(Xit,[S.IE.l],S.fluidviscosity));
        [S.IN,S.IE,S.SE] = solvehemodyn(S.IN,S.IE,S.SE);

        % Store
        mX_pulse(it,:)   = Xit';
        mQ_pulse(it,:)   = [S.IE.Q]';
        mWSS_pulse(it,:) = vout(calcshearstress(mQ_pulse(it,:), Xit, S.fluidviscosity));

        % Advance
        t_current = t_current + dt;
        it = it + 1;
    end

    %% --- 5. Add to combined plots ---
    if showRadius
        figure(fig_radius);
        mean_radius = mean(mX_pulse, 2);
        plot(tout_pulse, mean_radius, 'LineWidth', 1.5, ...
             'Color', colors_radius(ik,:), ...
             'DisplayName', ['k=' num2str(S.k)]);
    end

    if showFlow
        figure(fig_flow);
        mean_flow = mean(mQ_pulse, 2);
        plot(tout_pulse, mean_flow, 'LineWidth', 1.5, ...
             'Color', colors_flow(ik,:), ...
             'DisplayName', ['k=' num2str(S.k)]);
    end

    if showWSS
        figure(fig_wss);
        mean_wss = mean(mWSS_pulse, 2);
        plot(tout_pulse, mean_wss, 'LineWidth', 1.5, ...
             'Color', colors_wss(ik,:), ...
             'DisplayName', ['k=' num2str(S.k)]);
    end

    disp(['Finished simulation for k = ', num2str(S.k)]);

    %% --- 6. Cycle-Averaged WSS Analysis ---
    samples_per_cycle = round(pulse_period / dt);
    num_cycles = floor(length(tout_pulse) / samples_per_cycle);
    avgWSS = zeros(num_cycles, nElem);
    t_cycles = zeros(num_cycles, 1);

    for c = 1:num_cycles
        idx_start = (c-1)*samples_per_cycle + 1;
        idx_end = c*samples_per_cycle;
        avgWSS(c,:) = mean(mWSS_pulse(idx_start:idx_end,:), 1);
        t_cycles(c) = tout_pulse(idx_end);
    end

    meanWSS_all = mean(abs(avgWSS), 2) * 11;

    %% --- 7. Stabilization Detection ---
    deltaWSS = abs(diff(meanWSS_all)) ./ meanWSS_all(1:end-1);
    tolerance = 0.0001;
    disp(meanWSS_all);
    disp('diff');
    disp(diff(meanWSS_all));
    disp(deltaWSS);

    stable_idx = find(deltaWSS < tolerance, 1, 'first');
    if ~isempty(stable_idx)
        stable_time = t_cycles(stable_idx);
        fprintf('✅ Mean WSS stabilized for k = %.4f at t = %.2f s\n', S.k, stable_time);
    else
        fprintf('⚠️ WSS did not stabilize for k = %.4f within %.2f s\n', S.k, tend_pulse);
    end

    %% --- 8. Store for comparison & combined plot ---
    WSS_data(ik).k = S.k;
    WSS_data(ik).t_cycles = t_cycles;
    WSS_data(ik).meanWSS = meanWSS_all;
    %% --- Save WSS data for this k to CSV ---
    filename = sprintf('WSS_k_%.4f.csv', S.k);
    
    % Combine data into a table: time, mean WSS, and all vessel-wise WSS
    T = array2table([t_cycles, meanWSS_all, avgWSS], ...
        'VariableNames', [{'t_cycles', 'meanWSS'}, arrayfun(@(n) sprintf('Vessel_%d', n), 1:nElem, 'UniformOutput', false)]);
    
    % Write to CSV
    writetable(T, filename);
    fprintf('💾 Saved WSS data for k = %.4f to %s\n', S.k, filename);

    % Add to combined plot
    plot(t_cycles, meanWSS_all, 'LineWidth', 1.8, 'Color', colors(ik,:), ...
         'DisplayName', ['k = ', num2str(S.k)]);
end

%% === 9. Finalize combined plots ===
if showRadius
    figure(fig_radius);
    hold off;
    xlabel('Time [s]'); ylabel('Radius [m]');
    title('Radius under Pulsatile Pressure (all k values)');
    legend('show', 'Location', 'best');
    grid on;
end

if showFlow
    figure(fig_flow);
    hold off;
    xlabel('Time [s]'); ylabel('Flow [m^3/s]');
    title('Flow under Pulsatile Pressure (all k values)');
    legend('show', 'Location', 'best');
    grid on;
end

if showWSS
    figure(fig_wss);
    hold off;
    xlabel('Time [s]'); ylabel('Wall Shear Stress [Pa]');
    title('WSS under Pulsatile Pressure (all k values)');
    legend('show', 'Location', 'best');
    grid on;
end

%% === 10. Final combined WSS plot formatting ===
hold off;
xlabel('Time [s]');
ylabel('Mean WSS over cycle [Pa]');
title('Cycle-averaged WSS evolution for different k');
legend('show', 'Location', 'best');
grid on;

%% === 11. Save all results ===
save('WSS_data_allK.mat', 'WSS_data');
disp('💾 Saved all WSS data to WSS_data_allK.mat');