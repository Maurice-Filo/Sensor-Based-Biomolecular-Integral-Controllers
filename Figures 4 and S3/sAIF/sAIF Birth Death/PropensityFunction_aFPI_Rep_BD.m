function Prop = PropensityFunction_aFPI_Rep_BD(x, Parameters)
% Propensity for Network aFPI_Rep_BD
% 	 Species: 		 X = [X_1; Z_1; Z_2]
% 	 Reactions: 	R1:     X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		phi				-->     Z_1				[mu]
% 				    R3:		phi             -->     Z_2				[theta*X_1]
% 				    R4:		Z_1 + Z_2		--> 	0				[eta*Z_1*Z_2]
% 				    R5:		phi				--> 	X_1				[alpha / (1 + (Z_2/kappa)^n )]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
mu = Parameters.mu;
alpha = Parameters.alpha;
theta = Parameters.theta;
eta = Parameters.eta;
kappa = Parameters.kappa;
n = Parameters.n;

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
		alpha/(1+(Z_2/kappa)^n); ...
		];
end