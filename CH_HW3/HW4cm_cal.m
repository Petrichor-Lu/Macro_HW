function run_tau_loop()
    % parameters
    delta = 0.025;
    alpha = 0.3;
    beta  = 0.99;
    theta = 0.384;

    tau_grid = linspace(0, 1, 21);  % 21 grid-term
    n_grid = length(tau_grid);

    % result storage
    results = zeros(n_grid, 6); % c, l, k, w, r, T

    % initial guess
    x0 = [1, 0.3, 5, 1, 0.04];

    for i = 1:n_grid
        tau = tau_grid(i);

        % us fsolve solve
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10);
        x_ss = fsolve(@(x) ss_eqm(x, alpha, beta, delta, theta, tau), x0, options);

        % save results
        c = x_ss(1); l = x_ss(2); k = x_ss(3); w = x_ss(4); r = x_ss(5);
        T = tau * w * l;

        results(i, :) = [c, l, k, w, r, T];

        fprintf('\n[ tau = %.3f ]\n', tau);
        fprintf('Consumption (c)   = %.4f\n', c);
        fprintf('Labor (l)         = %.4f\n', l);
        fprintf('Capital (k)       = %.4f\n', k);
        fprintf('Wage (w)          = %.4f\n', w);
        fprintf('Interest rate (r) = %.4f\n', r);
        fprintf('Tax Revenue (T)   = %.4f\n', T);
    end

    % plot
    var_names = {'Consumption (c)', 'Labor (l)', 'Capital (k)', ...
                 'Wage (w)', 'Interest Rate (r)', 'Tax Revenue (T)'};

    for j = 1:6
        figure;
        plot(tau_grid, results(:, j), '-o', 'LineWidth', 2);
        xlabel('\tau (Labor Tax Rate)');
        ylabel(var_names{j});
        title(['\tau vs ', var_names{j}]);
        grid on;
    end
end

% function system
function F = ss_eqm(x, alpha, beta, delta, theta, tau)
    c = x(1);
    l = x(2);
    k = x(3);
    w = x(4);
    r = x(5);

    F(1) = (theta / c) * w * (1 - tau) - (1 - theta) / (1 - l);
    F(2) = beta * (1 + r - delta) - 1;
    F(3) = (1 - alpha) * k^alpha * l^(-alpha) - w;
    F(4) = alpha * k^(alpha - 1) * l^(1 - alpha) - r;
    F(5) = c + delta * k - k^alpha * l^(1 - alpha);
end
