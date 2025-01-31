function X_bar = aFPI_Deg_FixedPoint(Parameters, SupportingInput)
% u = alpha - gamma*z_2 * x_1/(x_1 + kappa_1)
%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
alpha = Parameters.alpha;
eta = Parameters.eta;
gamma = Parameters.gamma;
kappa_1 = Parameters.kappa_1;
r = mu / theta;

%% Compute Supporting Input
[u, x] = SupportingInput(Parameters, r);
T = x(1)/kappa_1 / (1 + x(1)/kappa_1);

%% Compute Fixed Point
z_1 = (mu*gamma*T) / (eta*(alpha-u));
z_2 = (alpha-u) / (gamma*T);

%% Stack Coordinates
X_bar = [x; z_1; z_2];
end

