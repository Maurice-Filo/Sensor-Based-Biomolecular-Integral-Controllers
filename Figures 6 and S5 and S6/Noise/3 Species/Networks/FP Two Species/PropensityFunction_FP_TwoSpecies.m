function Prop = PropensityFunction_FP_TwoSpecies(x, Parameters)
% Propensity Function for Two Species Process Controlled by
% Filtered Proportional Controller
% 	 Species: 		 X = [X_1; X_2; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
%                   R2:		 X_1				--> 	X_1 + X_2       [k_1*X_1]
%                   R3:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R4:		 0                  --> 	X_1				[alpha/(1 + Z_2/kappa)]
% 				    R5:		 X_2                --> 	X_2 + Z_2		[theta*X_2]
% 				    R6:		 Z_2                --> 	0				[delta*Z_2]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;
alpha = Parameters.alpha;
kappa = Parameters.kappa;
theta = Parameters.theta;
delta = Parameters.delta;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z_2 = x(3);

%% Propensities
Prop = [ ...
        gamma_1*X_1; ...
        k_1*X_1; ...
        gamma_2*X_2; ...
        alpha/(1 + Z_2/kappa); ...
        theta*X_2; ...
        delta*Z_2; ...
	];
end