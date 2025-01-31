function FP = ComputeFP_sAIF_BirthDeath(Parameters)
%% Extract Parameters
gamma_1 = Parameters.gamma_1;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
delta = Parameters.delta;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Solve Polynomial Equation
x_L = roots(flip([ ...
                    alpha^2*delta*eta*kappa^2; 
                    alpha*gamma_1*kappa*(delta^2 - 2*eta*kappa*delta + eta*mu); 
                   -gamma_1*kappa*(gamma_1*delta^2 - eta*gamma_1*kappa*delta + eta*gamma_1*mu + alpha*eta*theta); 
                   -gamma_1^2*theta*(delta - eta*kappa); ...
               ]));

%% Remove Negative and Complex Roots
x_L = x_L(x_L > 0 & imag(x_L) == 0);

%% Remove ROOTS that are not Feasible
x_L = x_L(x_L < alpha/gamma_1);

%% Check what is left
if length(x_L)~=1
    FP = NaN;
    return
end

%% Compute Remaining Coordinates of the Fixed Point
z_2 = kappa*(alpha/(gamma_1*x_L) - 1);
z_1 = (mu + delta*z_2 - theta*x_L)/delta;
FP = [x_L; z_1; z_2];

end