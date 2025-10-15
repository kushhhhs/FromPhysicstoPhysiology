%% === STEADY-STATE ADAPTATION ONLY ===
clear; clc; close all;

% === Choose k-values ===
k_values = [0.001, 0.0025, 0.005, 0.01];

% === Prepare figures ===
fig_mean = figure; hold on;
colors = lines(length(k_values));

% === Data storage structure ===
steady_data = struct();

for ik = 1:length(k_values)
    %% --- 1. Model setup ---
    S = struct();
    S.Assignment = 2;
    warning('off','MATLAB:singularMatrix')
    warning('off','MATLAB:nearlySingularMatrix')

    S = ModelPars(S);
    S = DefineTopology(S);
    [S.SE(find(S.sources)).sourceP] = deal(S.sourceP);
    [S.SE(find(~S.sources)).sourceP] = deal(S.sinkP);
    [S.IE.r] = vout([S.IE.r0]);
    [S.IE.G] = vout(conductance([S.IE.r],[S.IE.l],S.fluidviscosity));
    [S.IN,S.IE,S.SE] = solvehemodyn(S.IN,S.IE,S.SE);

    %% --- 2. Simulation parameters ---
    S.tend = 20000;                % total simulation time [s]
    S.k = k_values(ik);            % adaptation rate [1/s]
    S.WSSref = 5;                  % reference wall shear stress [Pa]
    options = odeset('RelTol',1e-4,'MaxStep',1);
    X0 = [S.IE.r0]';               % initial radius vector

    %% --- 3. Steady-state simulation ---
    [tout, mX] = ode45(@(t,X) adapt(t,X,S), [0 S.tend], X0, options);

    % Compute mean radius
    mean_r = mean(mX, 2);

    %% --- 4. Plot mean radius ---
    figure(fig_mean);
    plot(tout, mean_r, 'LineWidth', 1.5, ...
         'Color', colors(ik,:), ...
         'DisplayName', ['k = ' num2str(S.k)]);
    drawnow;

    %% --- 5. Store results ---
    steady_data(ik).k = S.k;
    steady_data(ik).time = tout;
    steady_data(ik).mean_radius = mean_r;
    steady_data(ik).final_radius = mX(end,:);

    %% --- 6. Print final radii in the command window ---
    fprintf('\n=== Final radii for k = %.4f ===\n', S.k);
    disp(mX(end,:));
    fprintf('Mean final radius: %.6e m\n\n', mean(mX(end,:)));

    disp(['✅ Finished steady-state run for k = ', num2str(S.k)]);
end

%% --- 7. Plot formatting ---
figure(fig_mean);
xlabel('Time [s]');
ylabel('Mean Radius [m]');
title('Steady-State Radius Adaptation for Different k Values');
legend('show', 'Location', 'best');
grid on;

%% --- 8. Save results ---
save('steady_radius_allK.mat', 'steady_data');
disp('💾 Saved steady-state adaptation data to steady_radius_allK.mat');
