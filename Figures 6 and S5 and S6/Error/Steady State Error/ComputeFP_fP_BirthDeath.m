function FP = ComputeFP_fP_BirthDeath(Parameters)
%% Extract Parameters
gamma_1 = Parameters.gamma_1;
theta = Parameters.theta;
delta = Parameters.delta;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Solve Polynomial Equation
x_L = (delta*kappa/theta/2) * (-1 + sqrt(1 + 4*alpha*theta/gamma_1/delta/kappa));

%% Compute Remaining Coordinates of the Fixed Point
z_2 = theta*x_L/delta;
FP = [x_L; z_2];

end