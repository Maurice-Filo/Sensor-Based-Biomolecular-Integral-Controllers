function Stability = ComputeStability_sAIF_BirthDeath(Parameters, FP)
%% Extract Parameters
gamma_1 = Parameters.gamma_1;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
delta = Parameters.delta;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Extract Coordinates
x_L = FP(1);
z_1 = FP(2);
z_2 = FP(3);

%% Compute Jacobian
A = [-gamma_1,      0,                  -alpha/(kappa*(z_2/kappa + 1)^2); ...
     0,            -delta - eta*z_2,    -eta*z_1; ...
     theta,        -eta*z_2,            -delta - eta*z_1];

%% Compute Eigenvalues
EIGENVALUES = eig(A);

%% Check Stability
Stability = ~any(real(EIGENVALUES) > 0);

end