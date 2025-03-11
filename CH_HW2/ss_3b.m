% PS#2 Question 3.(b) S.S.

% Define initial guess values
x0 = [0.01, 0.2, 500, 1, 2, 4, 3]; % [r, u, K, i, W, Y, C]

% Solve the nonlinear system of equations using fsolve
options = optimoptions('fsolve', 'Display', 'iter'); % Display iteration process
sol = fsolve(@steady_state_equations, x0, options);

% Store steady-state values in H matrix
H = [sol(7); sol(3); 1; sol(6); sol(4); sol(2); sol(5); sol(1)]; % [C; K; L; Y; I; U; W; R]

% Convert H to double precision to avoid compatibility issues
H = double(H);

H

% Save H in MATLAB v7 format to ensure Dynare compatibility
save('steady_state_values.mat', 'H', '-v7');

% Define the system of equations
function F = steady_state_equations(x)
    % Variables
    r = x(1);
    u = x(2);
    K = x(3);
    i = x(4);
    W = x(5);
    Y = x(6);
    C = x(7);

    % Parameters
    beta = 0.99;
    alpha = 0.3;
    delta_0 = 0;
    delta_1 = 0.03;
    omega = 2;

    % 7 steady-state equations
    F(1) = 0.99 * (r * u + 1 - 0.03 * u^2) - 1;  % Capital market equilibrium
    F(2) = r - 0.06 * u;                         % Return on capital
    F(3) = r - alpha * K^(-0.7) * u^(-0.7);      % Capital stock condition
    F(4) = i - (delta_0 + delta_1 * u^omega) * K; % Capital accumulation
    F(5) = W - (1 - alpha) * K^alpha * u^alpha;  % Wage equation
    F(6) = Y - K^alpha * u^alpha;                % Output equation
    F(7) = Y - C - i;                            % Resource constraint
end
