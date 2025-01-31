function Gains = aFPI_Rep_Params2Gains(Parameters, SupportingInput)
% u = alpha/(1+Z_2/kappa)

%% Extract Controller Parameters
theta = Parameters.theta;
eta = Parameters.eta;
kappa = Parameters.kappa;
alpha = Parameters.alpha;
mu = Parameters.mu;
r = mu/theta;
n = Parameters.n;

%% Fixed Point
X_bar = aFPI_Rep_FixedPoint(Parameters, SupportingInput);
Z_bar_1 = X_bar(2);
Z_bar_2 = X_bar(3);

[u_bar, ~] = SupportingInput(Parameters, r);

%% Compute Derivatives
sigma_1 = n*u_bar^2/alpha/kappa*(Z_bar_2/kappa)^(n-1);

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

