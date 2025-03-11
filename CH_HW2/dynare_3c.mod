// PS#2 Question 3.(c)

// 1. Declare endogenous and exogenous variables
var y c i k w r z u;
varexo e; 

// 2. Declare parameters
parameters beta alpha delta_0 delta_1 omega rho_z;

beta    = 0.99;   // Discount factor
alpha   = 0.3;    // Capital share
delta_0 = 0;      // Constant part of depreciation
delta_1 = 0.03;   // Depreciation coefficient
omega   = 2;      // Capital utilization parameter
rho_z   = 0.9;    // AR(1) coefficient for TFP shock

// 3. Define model equations
model;

// 1) Euler Equation
1 / c = beta * ((1 / c(+1)) * (1 - delta_0 + (omega - 1) * delta_1 * u(+1)^omega));

// 2) Capital Utilization
r = delta_1 * omega * u^(omega - 1);

// 3) Capital Accumulation
k = i + (1 - delta_0 - delta_1 * u^omega) * k(-1);

// 4) Resource Constraint
c + i = y;

// 5) Wage Equation
w = (1 - alpha) * exp(z) * (u * k(-1))^alpha;

// 6) Interest Rate Equation
r = alpha * exp(z) * (u * k(-1))^(alpha-1);

// 7) Production Function
y = exp(z) * (u * k(-1))^alpha;

// 8) TFP Shock Process
z = rho_z * z(-1) + e;

end;

// 4. Load steady-state values from MATLAB file
load steady_state_values.mat;

// 5. Set initial values
initval;
    z = 0;
    y = H(4);
    c = H(1);
    i = H(5);
    k = H(2);
    w = H(7);
    r = H(8);
    u = H(6);
end;

steady;

// 6. Define shocks
shocks;
    var e; stderr 0.01; // Standard deviation of TFP shock
end;

// 7. Compute solution and impulse response functions (IRFs)
stoch_simul(order=1, irf=20, periods=1000, hp_filter=1600) c i y k;
