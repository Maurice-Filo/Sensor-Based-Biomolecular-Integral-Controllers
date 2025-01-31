function X_bar = I_FixedPoint(Parameters, SupportingInput)
% u = k * z_1

%% Extract Controller Parameters
mu = Parameters.mu;
theta = Parameters.theta;
k = Parameters.k;
eta = Parameters.eta;
r = mu / theta;

%% Compute Supporting Input
[u, x] = SupportingInput(Parameters, r);

%% Compute Fixed Point
z_1 = u / k;
z_2 = mu / (eta*z_1);

%% Stack Coordinates
X_bar = [x; z_1; z_2];
end

