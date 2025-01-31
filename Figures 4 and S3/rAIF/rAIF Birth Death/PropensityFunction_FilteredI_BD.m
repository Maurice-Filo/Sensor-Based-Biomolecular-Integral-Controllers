function Prop = PropensityFunction_FilteredI_BD(x, Parameters)
% Propensity for Network FilteredI_BD
% 	 Species: 		 X = [X_1; Z_1; Z_2]
% 	 Reactions: 	R1:     X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		phi				-->     Z_1				[mu]
% 				    R3:		phi             -->     Z_2				[theta*X_1]
% 				    R4:		Z_1 + Z_2		--> 	0				[eta*Z_1*Z_2]
% 				    R5:		phi				--> 	X_1				[k*Z_1]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
k = Parameters.k;

%% Extract State Variables
X_1 = x(1);
Z_1 = x(2);
Z_2 = x(3);

%% Propensities
Prop = [ ...
		gamma_1*X_1; ...
		mu; ...
		theta*X_1; ...
		eta*Z_1*Z_2; ...
		k*Z_1; ...
		];
end