function FP = ComputeFP_rAIF_BirthDeath(Parameters)
%% Extract Parameters
gamma_1 = Parameters.gamma_1;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
delta = Parameters.delta;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Solve Polynomial Equation
x_L = roots(flip([-alpha^2*delta*mu; ...
                   alpha*gamma_1*(kappa*delta^2 + 2*mu*delta - eta*kappa*mu); ...
                  -delta^2*gamma_1^2*kappa + eta*delta*gamma_1^2*kappa^2 - mu*delta*gamma_1^2 + eta*mu*gamma_1^2*kappa + alpha*eta*theta*gamma_1*kappa; ...
                  -eta*gamma_1^2*kappa*theta]));

%% Remove Negative and Complex Roots
x_L = x_L(x_L > 0 & imag(x_L) == 0);

%% Remove ROOTS that are not Feasible
x_L = x_L(x_L < alpha/gamma_1);

%% Check what is left
if length(x_L)~=1
    FP = NaN;
    return;
end

%% Compute Remaining Coordinates of the Fixed Point
z_1 = (gamma_1*kappa*x_L)/(alpha - gamma_1*x_L);
z_2 = (delta*z_1 - mu + theta*x_L)/delta;
FP = [x_L; z_1; z_2];

end