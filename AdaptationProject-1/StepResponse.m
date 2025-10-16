%% === PARAMETER SWEEP: DIFFERENT ADAPTATION SPEEDS (k) ===
clear; clc; close all;

% === Plot display options ===
showRadius = true;
showFlow   = false;
showWSS    = false;

% === Define adaptation rate values ===
k_values = [0.001, 0.0025, 0.005, 0.01];

% === Preallocate storage ===
final_radii = cell(length(k_values),1);
data = struct();

% === Setup figures before loop ===
if showRadius
    fig_radius = figure; hold on; colors_radius = lines(length(k_values));
end
if showFlow
    fig_flow = figure; hold on; colors_flow = lines(length(k_values));
end
if showWSS
    fig_wss = figure; hold on; colors_wss = lines(length(k_values));
end

% === Combined WSS figure setup ===
figure; hold on; colors = lines(length(k_values));

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

    disp('💡 WSS just before step flow:');
    disp([S.IE.WSS]);

    %% --- 2. Simulation parameters ---
    S.tend = 400;                % steady adaptation duration [s]
    S.k = k_values(ik);          % adaptation rate [1/s]
    S.WSSref = 5;                % reference WSS [Pa]
    options = odeset('RelTol',1e-4,'MaxStep',1);
    X0 = [S.IE.r0]';

    %% --- 3. STEADY ADAPTATION RUN ---
    [tout, mX] = ode45(@(t,X) adapt(t,X,S), [0 S.tend], X0, options);
    final_radii{ik} = mX(end,:);

    %% --- 4. STEP RESPONSE SIMULATION ---
    X0_step = mX(end,:);
    t_current = 0; dt = 0.01;
    tend_step = 400;          % total time [s]
    step_time = 5;            % step occurs at t = 5 s
    step_amp = 0.2;           % 30% increase in inlet pressure

    tout_step = t_current:dt:tend_step;
    nElem = length(X0_step);
    mX_step   = zeros(length(tout_step), nElem);
    mQ_step   = zeros(length(tout_step), nElem);
    mWSS_step = zeros(length(tout_step), nElem);
    Xit = X0_step; it = 1;

    while t_current <= tend_step
        % Apply step input
        if t_current >= step_time
            stepP = S.sourceP * (1 + step_amp);
        else
            stepP = S.sourceP;
        end
        [S.SE(find(S.sources)).sourceP] = deal(stepP);
        [S.SE(find(~S.sources)).sourceP] = deal(S.sinkP);

        % Adaptation step
        tspan = [0 dt];
        [~, X_step] = ode45(@(t,X) adapt(t,X,S), tspan, Xit);
        Xit = X_step(end,:)';
        [S.IE.r] = vout(Xit);
        [S.IE.G] = vout(conductance(Xit,[S.IE.l],S.fluidviscosity));
        [S.IN,S.IE,S.SE] = solvehemodyn(S.IN,S.IE,S.SE);

        % Store results
        mX_step(it,:)   = Xit';
        mQ_step(it,:)   = [S.IE.Q]';
        mWSS_step(it,:) = vout(calcshearstress(mQ_step(it,:), Xit, S.fluidviscosity));

        t_current = t_current + dt;
        it = it + 1;
    end

    %% --- 5. Plotting (if enabled) ---
    if showRadius
        figure(fig_radius);
        plot(tout_step, mean(mX_step,2), 'LineWidth', 1.5, ...
             'Color', colors_radius(ik,:), ...
             'DisplayName', ['k = ', num2str(S.k)]);
    end
    if showFlow
        figure(fig_flow);
        plot(tout_step, mean(mQ_step,2), 'LineWidth', 1.5, ...
             'Color', colors_flow(ik,:), ...
             'DisplayName', ['k = ', num2str(S.k)]);
    end
    if showWSS
        figure(fig_wss);
        plot(tout_step, mean(mWSS_step,2), 'LineWidth', 1.5, ...
             'Color', colors_wss(ik,:), ...
             'DisplayName', ['k = ', num2str(S.k)]);
    end

    disp(['✅ Finished step simulation for k = ', num2str(S.k)]);

    %% --- 6. Stabilization detection ---
    meanWSS_all = mean(abs(mWSS_step), 2);
    deltaWSS = abs(diff(meanWSS_all)) ./ meanWSS_all(1:end-1);
    tolerance = 0.01;  % 1% threshold

    stable_idx = find(deltaWSS < tolerance, 1, 'first');
    if ~isempty(stable_idx)
        stable_time = tout_step(stable_idx);
        fprintf('✅ Mean WSS stabilized for k = %.4f at t = %.2f s\n', S.k, stable_time);
    else
        fprintf('⚠️ WSS did not stabilize for k = %.4f within %.2f s\n', S.k, tend_step);
    end

    %% --- 7. Save to structure and CSV ---
    data(ik).k = S.k;
    data(ik).tout = tout_step;
    data(ik).meanWSS = meanWSS_all;
    data(ik).Flow = mQ_step;
    data(ik).Radius = mX_step;
    filename = sprintf('StepResponse_WSS_k_%.4f.csv', S.k);
    T = array2table([tout_step', meanWSS_all, mWSS_step], ...
        'VariableNames', [{'time_s','meanWSS'}, ...
         arrayfun(@(n) sprintf('Vessel_%d', n), 1:nElem, 'UniformOutput', false)]);
    writetable(T, filename);
    fprintf('💾 Saved step WSS data for k = %.4f to %s\n', S.k, filename);
    filename_Flow = sprintf('StepResponse_Flow_k_%.4f.csv', S.k);
    T_Flow = array2table([tout_step', mean(mQ_step,2), mQ_step], ...
        'VariableNames', [{'time_s','meanFlow'}, ...
         arrayfun(@(n) sprintf('Vessel_%d', n), 1:nElem, 'UniformOutput', false)]);
    writetable(T_Flow, filename_Flow);
    fprintf('💾 Saved Flow data for k = %.4f to %s\n', S.k, filename_Flow);
    filename_Radius = sprintf('StepResponse_Radius_k_%.4f.csv', S.k);
    T_Radius = array2table([tout_step', mean(mX_step,2), mX_step], ...
        'VariableNames', [{'time_s','meanRadius'}, ...
         arrayfun(@(n) sprintf('Vessel_%d', n), 1:nElem, 'UniformOutput', false)]);
    writetable(T_Radius, filename_Radius);
    fprintf('💾 Saved Radius data for k = %.4f to %s\n', S.k, filename_Radius);
    % % Combined plot
    % plot(tout_step, meanWSS_all, 'LineWidth', 1.8, ...
    %     'Color', colors(ik,:), 'DisplayName', ['k = ', num2str(S.k)]);
end

%% --- 8. Finalize plots ---
if showRadius
    figure(fig_radius); hold off;
    xlabel('Time [s]'); ylabel('Radius [m]');
    title('Radius evolution under step flow');
    legend('show'); grid on;
end
if showFlow
    figure(fig_flow); hold off;
    xlabel('Time [s]'); ylabel('Flow [m^3/s]');
    title('Flow evolution under step flow');
    legend('show'); grid on;
end
if showWSS
    figure(fig_wss); hold off;
    xlabel('Time [s]'); ylabel('WSS [Pa]');
    title('Wall shear stress under step flow');
    legend('show'); grid on;
end

xlabel('Time [s]');
ylabel('Mean WSS [Pa]');
title('Mean WSS evolution for different adaptation rates (Step Flow)');
legend('show'); grid on;

save('StepResponse_WSS_allK.mat', 'data');
disp('💾 Saved all WSS data to StepResponse_WSS_allK.mat');
