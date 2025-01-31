function Gains = aFPI_Deg_Params2Gains(Parameters, SupportingInput)
% u = alpha - gamma*z_2 * x_1/(x_1 + kappa_1)

%% Extract Controller Parameters
% mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
kappa_1 = Parameters.kappa_1;
gamma = Parameters.gamma;
% r = mu / theta;

%% Fixed Point
X_bar = aFPI_Deg_FixedPoint(Parameters, SupportingInput);
X_bar_1 = X_bar(1);
Z_bar_1 = X_bar(2);
Z_bar_2 = X_bar(3);

T = (X_bar_1/kappa_1) / (1+X_bar_1/kappa_1);

%% Compute Derivatives
sigma_1 = gamma * T;

%% Compute cut-off frequency
omega_c = eta*(Z_bar_1+Z_bar_2);

%% Compute PID Gains
K_P = sigma_1/omega_c;
K_I = (sigma_1*Z_bar_2)/(Z_bar_1+Z_bar_2);
K_S = theta;

%% Construct Gains Structure
Gains.K_P = K_P;
Gains.K_I = K_I;
Gains.omega_c = omega_c;
Gains.K_S = K_S;

end

