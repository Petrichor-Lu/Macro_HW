%% Question (e): Effect of varying kappa on equilibrium consumption and prices
clear               
clc                
format compact     

% --- Model parameters ---
omega = 0.3;        % CES weight on tradable goods
eta = 0.205;        % CES elasticity parameter
b0 = -1.0064;       % Initial bond holdings
yN = 1;             % Nontradable endowment (fixed)
yTt = 1;            % Steady-state tradable endowment
shock = 0.05;       % Fixed 5% negative shock to tradable income
yT0 = (1 - shock) * yTt;

% --- Convergence controls ---
tol = 1e-8;         % Tolerance for fixed-point convergence
max_iter = 1000;    % Maximum iterations allowed

% --- Define range of kappa values ---
kappa_grid = 0.1:0.01:0.3;
n_k = length(kappa_grid);

% --- Store results for each kappa ---
cT_kappa = zeros(1, n_k);   % Equilibrium c_0^T
pN_kappa = zeros(1, n_k);   % Equilibrium p_0^N

% --- Loop over each value of kappa ---
for i = 1:n_k
    kappa = kappa_grid(i);

    % Initial guess for c0 based on MM curve with old price
    pNt = 2.2206;  % Use steady-state price as initial reference
    c0_old = yT0 * (1 + kappa) + kappa * pNt + b0;
    
    diff = 1;
    iter = 0;

    % Fixed-point iteration: alternate between BB and MM curve
    while diff > tol && iter < max_iter
        % Update p_0^N from BB curve
        pN0_new = ((1 - omega) / omega) * c0_old^(1 + eta);
        
        % Update c_0^T from MM curve
        c0_new = yT0 * (1 + kappa) + kappa * pN0_new + b0;
        
        % Check convergence
        diff = abs(c0_new - c0_old);
        c0_old = c0_new;
        iter = iter + 1;
    end

    % Store results
    cT_kappa(i) = c0_new;
    pN_kappa(i) = pN0_new;
end

% --- Plot equilibrium values against kappa ---
figure;
plot(kappa_grid, cT_kappa, '-*', 'LineWidth', 2, ...
     'Color', [0 0.4470 0.7410], 'MarkerSize', 8, ...
     'DisplayName', 'c_0^T'); hold on;

plot(kappa_grid, pN_kappa, '-d', 'LineWidth', 2, ...
     'Color', [0.5 0 0], 'MarkerSize', 8, ...
     'DisplayName', 'p_0^N');


xlabel('\kappa (Borrowing limit parameter)');
ylabel('Equilibrium values');
title('Effect of \kappa on Consumption and Nontradable Prices');
legend('Location','northeast');
grid on;
