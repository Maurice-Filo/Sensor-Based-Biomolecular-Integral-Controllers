function X_bar = aFPI_Rep_FixedPoint(Parameters, SupportingInput)
% u = alpha/(1+(Z_2/kappa)^n)
%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
alpha = Parameters.alpha;
eta = Parameters.eta;
kappa = Parameters.kappa;
n = Parameters.n;
r = mu / theta;

%% Compute Supporting Input
[u_bar, x] = SupportingInput(Parameters, r);

%% Compute Fixed Point
z_1 = mu/(eta*kappa*(alpha/u_bar-1)^(1/n));
z_2 = kappa*(alpha/u_bar-1)^(1/n);

%% Stack Coordinates
X_bar = [x; z_1; z_2];
end

