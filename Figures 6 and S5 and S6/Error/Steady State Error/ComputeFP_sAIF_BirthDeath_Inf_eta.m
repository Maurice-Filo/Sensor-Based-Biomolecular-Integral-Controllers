function x_L = ComputeFP_sAIF_BirthDeath_Inf_eta(Parameters)
%% Extract Parameters
gamma_1 = Parameters.gamma_1;
mu = Parameters.mu;
theta = Parameters.theta;
delta = Parameters.delta;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Output
x_L = min(alpha/gamma_1, (mu - delta*kappa + sqrt((mu - delta*kappa)^2 + 4*theta*delta*kappa*alpha/gamma_1) ) / (2*theta));

end