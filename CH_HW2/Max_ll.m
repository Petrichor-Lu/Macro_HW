%% **Problem 1: MLE estimation using Kalman Filter**
clear
clc

% **1. Load data and demean y**
data = readtable('PS1.csv');  % Load data file
y = data.y;                   % Extract column y
y_tilde = y - mean(y);        % Compute demeaned y (subtract mean)

% **2. Set initial parameter values**
% X represents [rho, sigma_sq], where:
% - rho: State transition parameter
% - sigma_sq: Process noise variance
X_init = [0, 1];   % Initial guess: rho = 0, sigma_sq = 1

% **3. Maximize the log-likelihood using fminsearch**
% - fminsearch automatically adjusts X to minimize -kalman_ll(y_tilde, X)
% - Since fminsearch minimizes a function, we negate kalman_ll
mle_est = fminsearch(@(X) -kalman_ll(y_tilde, X), X_init);
disp('Maximum Likelihood Estimation (MLE) results:')
disp(mle_est)

% **4. Compute OLS estimation**
% OLS estimates a simple autoregressive model: y_t = rho * y_{t-1} + ε_t
y_dep = y_tilde(2:end);    % y_t (starting from t=2)
y_lag = y_tilde(1:end-1);  % y_{t-1} (starting from t=1)

% OLS estimation of rho (Least Squares)
rho_hat = (y_lag' * y_lag) \ (y_lag' * y_dep);

% Compute residuals
resid = y_dep - rho_hat * y_lag;
% Compute error variance sigma^2
sigma2_hat = (resid' * resid) / (length(resid) - 1);

% Output OLS estimation
ols_est = [rho_hat, sigma2_hat];
disp('Ordinary Least Squares (OLS) results:')
disp(ols_est)
