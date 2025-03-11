%% **Kalman Filter Log-Likelihood Calculation Function**
function LL = kalman_ll(y, X)
    % **1. Parse parameters**
    rho = X(1);       % State transition parameter (measuring the relationship between y_t and y_{t-1})
    sigma_sq = X(2);  % Process noise variance (measuring state uncertainty)

    % **2. Define the state-space model**
    % This is the core state-space representation of the Kalman filter:
    % - a_t = T * a_{t-1} + ϵ_t  (Transient Equation: State transition equation)
    % - y_t = Z * a_t + ε_t      (Measurement Equation: Observation equation)
    
    T = rho;        % State transition matrix (T affects how the state evolves)
    R = 1;          % Selection matrix (R affects how the state is influenced by noise)
    Q = sigma_sq;   % Process noise variance (affects the random fluctuations of the state)
    Z = 1;          % Measurement matrix (determines how the observation variable y_t depends on the state variable)
    H = 0;          % Measurement noise variance (assuming no measurement noise)

    N = length(y);  % Number of observations

    % **3. Initialize Kalman filter variables**
    a = 0;          % Initial state estimate (usually set to 0)
    P = eye(1);     % Initial state variance (uncertainty, set as an identity matrix)

    % **4. Compute predictions for the first time step t=1**
    % This is the first step of Kalman recursion, computing the first time step's prediction and uncertainty:
    a_pred = T * a;                           
    P_pred = T * P * T' + R * Q * R';

    % **Compute predicted observation value**
    y_pred = Z * a_pred;                       
    F = Z * P_pred * Z' + H;  % Prediction error variance (measuring prediction reliability)

    % **Compute prediction error (residual)**
    v = y(1) - y_pred;  

    % **Compute the log-likelihood for the first time step**
    LL = -0.5 * ( log(2 * pi) + log(F) + (v^2) / F );

    % **5. Compute Kalman gain and update state estimates**
    K = (P_pred * Z') / F;  
    a = a_pred + K * v;  
    P = P_pred - K * Z * P_pred;

    % **6. Recursive calculation for all time steps t=2 to t=N**
    for t = 2:N
        % **(a) Predict the state variable**
        % Predict the current state based on the previous updated state a:
        a_pred = T * a;                        
        P_pred = T * P * T' + R * Q * R';

        % **(b) Compute the predicted observation value**
        y_pred = Z * a_pred;  % Predict y_t
        F = Z * P_pred * Z' + H;  % Compute prediction error variance F_t

        % **(c) Compute the prediction error (residual)**
        v = y(t) - y_pred;  % Compute the prediction error for the current time step

        % **(d) Accumulate log-likelihood**
        % Here, `LL` accumulates the log-likelihood values across all time steps t:
        LL = LL - 0.5 * ( log(2*pi) + log(F) + (v^2)/F );

        % **(e) Compute Kalman gain**
        K = (P_pred * Z') / F;  

        % **(f) Update state estimates**
        a = a_pred + K * v;  
        P = P_pred - K * Z * P_pred;
    end
end
