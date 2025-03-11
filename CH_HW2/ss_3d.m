% PS#2 Question 3.(c) S.S. (Updated with original code style)

% Define parameter values
beta    = 0.99;
alpha   = 0.3;
delta_0 = 0.025;  % Fixed depreciation rate
rho_z   = 0.9;    % AR(1) coefficient for productivity shock
l       = 1;      % Inelastic labor supply

% Define initial guess values
x0 = [0.01, 500, 1, 2, 4, 3, 0]; % [r, K, Y, W, I, C, Z]

% Solve the nonlinear system of equations using fsolve
options = optimoptions('fsolve', 'Display', 'iter'); % Display iteration process
sol = fsolve(@steady_state_equations, x0, options);

% Store steady-state values in H matrix
H = [sol(6); sol(2); l; sol(3); sol(5); sol(4); sol(1); sol(7)]; % [C; K; L; Y; I; W; R; Z]

% Convert H to double precision to avoid compatibility issues
H = double(H);

% Display steady-state values
fprintf('Steady-State Solution:\n');
fprintf('r = %.6f\n', sol(1));
fprintf('K = %.6f\n', sol(2));
fprintf('Y = %.6f\n', sol(3));
fprintf('W = %.6f\n', sol(4));
fprintf('I = %.6f\n', sol(5));
fprintf('C = %.6f\n', sol(6));
fprintf('Z = %.6f\n', sol(7));

% Save H in MATLAB v7 format to ensure Dynare compatibility
save('steady_state_values_3d.mat', 'H', '-v7');

% Define the system of equations
function F = steady_state_equations(x)
    % Variables
    r = x(1);
    K = x(2);
    Y = x(3);
    W = x(4);
    I = x(5);
    C = x(6);
    Z = x(7); % TFP shock

    % Parameters
    beta = 0.99;
    alpha = 0.3;
    delta_0 = 0.025;
    rho_z = 0.9;
    l = 1;  % Inelastic labor supply

    % 7 steady-state equations
    F(1) = r - (delta_0 + 0.01); % Net capital return assumption
    F(2) = K - (r / alpha)^(1 / (alpha - 1)); % Capital return condition
    F(3) = Y - exp(Z) * K^alpha; % Output equation with productivity shock
    F(4) = W - (1 - alpha) * exp(Z) * K^alpha; % Wage equation
    F(5) = I - delta_0 * K; % Capital accumulation
    F(6) = C - (Y - I); % Resource constraint
    F(7) = Z - rho_z * Z; % Steady-state assumption for Z (Z = 0 in steady state)
end