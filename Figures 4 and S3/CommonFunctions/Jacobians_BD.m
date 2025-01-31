function [A, B, W, C] = Jacobians_BD(Parameters, r)
%% Extract Plant Parameters
gamma_1 = Parameters.gamma_1;
k_1 = Parameters.k_1;

%% Compute Fixed Point
[~, x] = SupportingInput_BD(Parameters, r);
X_bar_1 = x(1);

%% Compute Plant Jacobians
A = [-gamma_1, 0; k_1, -gamma_2];
B = [1; 0];
W = [0; X_bar_1];
C = [0, 1];
end

