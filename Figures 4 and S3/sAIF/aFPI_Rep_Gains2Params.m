function [Parameters, Feasible] = aFPI_Rep_Gains2Params(Gains, Parameters, r, SupportingInput)
%% Extract Gains
Feasible = true;
K_P = Gains.K_P;
K_I = Gains.K_I;
omega_0 = Gains.omega_c;

%% Extract Controller Parameters
mu = Parameters.mu;
n = Parameters.n;

%% Compute Supporting Input
[u_bar, ~] = SupportingInput(Parameters, r);

%% Compute kappa
kappa = mu/(omega_0-K_I/K_P)*(n*u_bar/mu*(omega_0-K_I/K_P)/omega_0/K_P-1)^(1/n);

%% Compute alpha
alpha = u_bar/(1-mu/n/u_bar*omega_0*K_P/(omega_0-K_I/K_P));

%% Compute eta
eta = ((K_I/K_P)*(omega_0-(K_I/K_P))) / mu;

%% Store Computed Parameters
Parameters.kappa = kappa;
Parameters.alpha = alpha;
Parameters.eta = eta;
end

