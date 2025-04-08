% Initial guess [c, l, k, w, r, tau_k]
x0 = [0.6, 0.3, 6.5, 1.7, 0.035, 0.4];

% Solve
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun',1e-10, 'TolX',1e-10);
[x_ss, fval, exitflag] = fsolve(@ss_capital_tax, x0, options);

% Display results
c     = x_ss(1);
l     = x_ss(2);
k     = x_ss(3);
w     = x_ss(4);
r     = x_ss(5);
tau_k = x_ss(6);

fprintf('Steady state under capital tax:\n');
fprintf('c      = %.4f\n', c);
fprintf('l      = %.4f\n', l);
fprintf('k      = %.4f\n', k);
fprintf('w      = %.4f\n', w);
fprintf('r      = %.4f\n', r);
fprintf('tau_k  = %.4f\n', tau_k);

% Welfare
theta = 0.384;
beta  = 0.99;
W = theta * log(c) + (1 - theta) * log(1 - l);
fprintf('Welfare = %.6f\n', W);
