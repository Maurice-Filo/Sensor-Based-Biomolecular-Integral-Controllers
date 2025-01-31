function [Parameters, Feasible] = I_Gains2Params(Gains, Parameters, r, SupportingInput)
%% Extract Gains
K_I = Gains.K_I;
omega_0 = Gains.omega_0;

%% Extract Parameters
mu = Parameters.mu;

%% Supporting Input
[u_bar, ~] = SupportingInput(Parameters, r);

%% Compute eta
eta = omega_0 * K_I / u_bar;

%% Compute k
k = ((u_bar * omega_0) / (2 * mu)) * (1 - sqrt(1 - (4*mu/u_bar) * (K_I/omega_0)) );

%% Store Computed Parameters
Parameters.k = k;
Parameters.eta = eta;

%% Feasibility
if omega_0 > 4 * mu * K_I / u_bar
    Feasible = 1;
else
    Feasible = 0;
end

end

