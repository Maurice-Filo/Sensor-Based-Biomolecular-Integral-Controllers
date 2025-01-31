function Prop = PropensityFunction_FP_BirthDeath(x, Parameters)
% Propensity Function for Birth-Death Process Controlled by
% Filtered Proportional Controller
% 	 Species: 		 X = [X_1; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		 0                  --> 	X_1				[alpha/(1 + Z_2/kappa)]
% 				    R3:		 X_1                --> 	X_1 + Z_2		[theta*X_1]
% 				    R4:		 Z_2                --> 	0				[delta*Z_2]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
alpha = Parameters.alpha;
kappa = Parameters.kappa;
theta = Parameters.theta;
delta = Parameters.delta;

%% Extract State Variables
X_1 = x(1);
Z_2 = x(2);

%% Propensities
Prop = [ ...
        gamma_1*X_1; ...
        alpha/(1 + Z_2/kappa); ...
        theta*X_1; ...
        delta*Z_2; ...
	];
end