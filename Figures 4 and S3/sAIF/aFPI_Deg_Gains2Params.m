function [Parameters, Feasible] = aFPI_Deg_Gains2Params(Gains, Parameters, r, SupportingInput)
%% Extract Gains
Feasible = true;
K_P = Gains.K_P;
K_I = Gains.K_I;
omega_c = Gains.omega_c;
K_S = Gains.K_S;

%% Extract Controller Parameters
mu = Parameters.mu;
kappa_1 = Parameters.kappa_1;

%% Compute Supporting Input
[u_bar, x] = SupportingInput(Parameters, r);
T = (x(1)/kappa_1) / (1 + x(1)/kappa_1);

%% Compute gamma
gamma = K_P * omega_c / T;

%% Compute alpha
alpha = u_bar + (mu*omega_c*K_P) / (omega_c-K_I/K_P);

%% Compute eta
eta = ((K_I/K_P)*(omega_c-(K_I/K_P))) / mu;

%% Store Computed Parameters
Parameters.gamma = gamma;
Parameters.alpha = alpha;
Parameters.eta = eta;
end

