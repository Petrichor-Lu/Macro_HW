// PS#2 Question 3.(c) Dynare Model 

// Declare endogenous and exogenous variables
var y c i k w r z;
varexo e;  // Exogenous technology shock

// Define parameters
parameters beta alpha delta_0 rho_z;

beta    = 0.99;     // Discount factor
alpha   = 0.3;      // Capital share in production
delta_0 = 0.025;    // Fixed depreciation rate
rho_z   = 0.9;      // AR(1) coefficient of TFP shock

// Define model equations
model;

    // Euler equation
    1 / c = beta * ((1 / c(+1)) * (1 - delta_0 + r(+1)));

    // Capital accumulation equation
    k = i + (1 - delta_0) * k(-1);

    // Resource constraint
    c + i = y;

    // Wage equation
    w = (1 - alpha) * exp(z) * k(-1)^alpha;

    // Interest rate equation
    r = alpha * exp(z) * k(-1)^(alpha - 1);

    // Production function
    y = exp(z) * k(-1)^alpha;

    // Shock process
    z = rho_z * z(-1) + e;

end;

// Load steady-state values from MATLAB
run ss_3d.m 

// Assign steady-state values to Dynare variables
initval;
    z = H(8,1);
    y = H(4,1);
    c = H(1,1);
    i = H(5,1);
    k = H(2,1);
    w = H(6,1);
    r = H(7,1);
end;

steady;

// Define shock process
shocks;
    var e; stderr 0.01; // Standard deviation of the TFP shock
end;

// Compute impulse response functions
stoch_simul(order=1, irf=20, periods=1000, hp_filter=1600) c i y k;
