%% === Parameters ===
clear; clc;

% Model parameters
alpha = 0.3;
delta = 0.025;
beta  = 0.99;
theta = 0.384;

tau_old = 0.1;
tau_new = 0.2;

% Known steady states (from HW3 output)
k0 = 7.1400;
k1 = 6.5906;
l0 = 0.3331;
l1 = 0.3075;
c0 = 0.6569;
c1 = 0.6064;
y0 = k0^alpha * l0^(1 - alpha);
y1 = k1^alpha * l1^(1 - alpha);

% Time and iteration settings
T = 70;
lambda = 0.5;
tol = 1e-8;
max_iter = 500;

%% === Step 1: Initial guess k_hat ===
k_hat = zeros(T+2,1);
k_hat(1) = k0;
for t = 2:T
    k_hat(t) = 0.9 * k_hat(t-1) + 0.1 * k1;
end
k_hat(T+1) = k1;
k_hat(T+2) = k1;

%% === Step 2: Labor FOC solver ===
solve_l = @(kt, kt1) fzero(@(lt) ...
    ((1 - theta)/(1 - lt)) - ...
    ((theta / (kt^alpha * lt^(1 - alpha) - kt1 + (1 - delta)*kt)) ...
    * (1 - tau_new) * (1 - alpha) * kt^alpha * lt^(-alpha)), [1e-3, 0.999]);

%% === Step 3: Time Path Iteration ===
for iter = 1:max_iter
    k_tilde = update_k(k_hat, T, k0, k1, alpha, delta, beta, theta, solve_l);
    diff = norm(k_tilde(1:T) - k_hat(1:T));

    if diff < tol
        break;
    end
    k_hat(1:T) = lambda * k_tilde(1:T) + (1 - lambda) * k_hat(1:T);
end

fprintf("Converged in %d iterations.\n", iter);
k_path = k_tilde(1:T+1);

%% === Step 4: Recover variables (l, c, y, w, r, T) ===
l_path = zeros(T+1,1);
c_path = zeros(T+1,1);
y_path = zeros(T+1,1);
w_path = zeros(T+1,1);
r_path = zeros(T+1,1);
tax_path = zeros(T+1,1);

% === Time 1: Use steady state labor and tau_old ===
l_path(1) = l0;
k_now = k_path(1);
k_next = k_hat(2);
y_path(1) = k_now^alpha * l0^(1 - alpha);
c_path(1) = y_path(1) - k_next + (1 - delta)*k_now;
w_path(1) = (1 - alpha) * k_now^alpha * l0^(-alpha);
r_path(1) = alpha * k_now^(alpha - 1) * l0^(1 - alpha);
tax_path(1) = tau_old * w_path(1) * l0;

% === Time 2~T+1: solve labor under tau_new ===
for t = 2:T+1
    try
        lt = solve_l(k_path(t), k_hat(t+1));
    catch
        warning('solve_l failed at t = %d. Using fallback l = %.3f', t, l1);
        lt = l1;
    end
    kt = k_path(t);
    kt1 = k_hat(t+1);
    yt = kt^alpha * lt^(1 - alpha);
    ct = yt - kt1 + (1 - delta) * kt;
    wt = (1 - alpha) * kt^alpha * lt^(-alpha);
    rt = alpha * kt^(alpha - 1) * lt^(1 - alpha);
    Tt = tau_new * wt * lt;

    l_path(t) = lt;
    y_path(t) = yt;
    c_path(t) = ct;
    w_path(t) = wt;
    r_path(t) = rt;
    tax_path(t) = Tt;
end

%% === Step 5: Plot results ===
t = 0:T;
plot_var(t, k_path, 'Capital (k_t)');
plot_var(t, l_path, 'Labor (l_t)');
plot_var(t, c_path, 'Consumption (c_t)');
plot_var(t, y_path, 'Output (y_t)');
plot_var(t, w_path, 'Wage (w_t)');
plot_var(t, r_path, 'Return on Capital (r_t)');
plot_var(t, tax_path, 'Tax Revenue (T_t)');

%% === Function: update_k ===
function k_tilde = update_k(k_hat, T, k0, k1, alpha, delta, beta, theta, solve_l)
    k_tilde = zeros(T+2,1);
    k_tilde(1) = k0;
    k_tilde(T+2) = k1;

    for t = 1:T
        h = @(kt1) ...
            (theta / (k_hat(t)^alpha * solve_l(k_hat(t), kt1)^(1 - alpha) ...
            - kt1 + (1 - delta)*k_hat(t))) ...
            - beta * theta / (kt1^alpha * solve_l(kt1, k_hat(t+2))^(1 - alpha) ...
            - k_hat(t+2) + (1 - delta)*kt1) ...
            * (1 + alpha * kt1^(alpha - 1) * solve_l(kt1, k_hat(t+2))^(1 - alpha) - delta);

        try
            k_tilde(t+1) = fzero(h, [0.1 * k1, 5 * k1]);
        catch
            warning('fzero failed at t = %d. Using fallback guess.', t);
            k_tilde(t+1) = 0.9 * k_hat(t) + 0.1 * k1;
        end
    end
end

%% === Function: Plotting helper ===
function plot_var(t, series, labelname)
    figure;
    plot(t, series, 'LineWidth', 2);
    title(labelname);
    xlabel('Time'); ylabel(labelname);
    grid on;
end
