function [u, x] = SupportingInput_BD(Parameters, r)
%% Extract Plant Parameters
gamma_1 = Parameters.gamma_1;

%% Compute Plant Fixed Point
x = zeros(1,1);
x(1) = r;

%% Compute Supporting Input
u = gamma_1 * x(1);

end

