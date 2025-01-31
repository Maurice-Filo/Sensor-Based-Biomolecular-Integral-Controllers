function [u, x] = SupportingInput_GeneExp(Parameters, r)
%% Extract Plant Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
beta = Parameters.beta;

%% Compute Plant Fixed Point
x = zeros(2,1);
x(2) = r;
x(1) = gamma_2 * x(2) / beta;

%% Compute Supporting Input
u = gamma_1 * x(1);

end

