function Gains = I_Params2Gains(Parameters, SupportingInput)
% u = k * z_1

%% Extract Controller Parameters
k = Parameters.k;
theta = Parameters.theta;
eta = Parameters.eta;

%% Compute Partial Derivatives of Actuation
sigma_1 = k;

%% Fixed Point
X_bar = I_FixedPoint(Parameters, SupportingInput);
Z_bar_1 = X_bar(2);
Z_bar_2 = X_bar(3);

%% Compute PID Gains
K_I = sigma_1 * Z_bar_1 / (Z_bar_1 + Z_bar_2);
K_S = theta;
omega_0 = eta * (Z_bar_1 + Z_bar_2);
K_F = sigma_1 / omega_0;

%% Construct Gains Structure
Gains.K_I = K_I;
Gains.K_S = K_S;
Gains.omega_0 = omega_0;
Gains.K_F = K_F;

end

