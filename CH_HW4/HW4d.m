%% Problem (d): Amplification under Different Shock Sizes

clear               
clc               
format compact      

% --- Parameter setup ---
omega = 0.3;
kappa = 0.3;
eta = 0.205;
b0 = -1.0064;
yN = 1;
pNt = 2.2206;   % From steady state in part (b)
yTt = 1;        % Original tradable endowment

% --- Convergence controls ---
tol = 1e-8;     % Convergence threshold
max_iter = 1000; % Maximum number of iterations

% --- Shock grid: from 5% to 20% ---
shock_grid = 0.05:0.01:0.20;
n_shock = length(shock_grid);

% --- Initialize result vectors ---
cT_path = zeros(1, n_shock);   % Store equilibrium c_0^T for each shock
pN_path = zeros(1, n_shock);   % Store equilibrium p_0^N for each shock

% --- Main loop: solve for equilibrium under each shock level ---
for i = 1:n_shock
    
    % Tradable endowment under current shock
    s = shock_grid(i);
    yT0 = (1 - s) * yTt;

    % Initial guess for c0 from MM curve
    c0_old = yT0 * (1 + kappa) + kappa * pNt + b0;
    diff = 1;
    iter = 0;

    % Fixed-point iteration: alternate updates of pN0 and c0 until convergence
    while diff > tol && iter < max_iter
        % Update pN0: from BB curve
        pN0_new = ((1 - omega) / omega) * c0_old^(1 + eta);
        
        % Update c0: from MM curve (resource + borrowing constraint)
        c0_new = yT0 * (1 + kappa) + kappa * pN0_new + b0;
        
        % Check convergence
        diff = abs(c0_new - c0_old);
        c0_old = c0_new;
        iter = iter + 1;
    end

    % Store equilibrium results
    cT_path(i) = c0_new;
    pN_path(i) = pN0_new;
end

% --- Plotting the result ---
figure;
plot(shock_grid, cT_path, '-*', 'LineWidth', 2, ...
     'Color', [0 0.4470 0.7410], 'MarkerSize', 8, ...
     'DisplayName', 'c_0^T'); hold on;

plot(shock_grid, pN_path, '-d', 'LineWidth', 2, ...
     'Color', [0.5 0 0], 'MarkerSize', 8, ...
     'DisplayName', 'p_0^N');

xlabel('s: Tradable endowment shock');
ylabel('Equilibrium values');
legend({'c_0^T','p_0^N'}, 'Location','northeast');
title('Amplification under Increasing Shocks');
grid on;
