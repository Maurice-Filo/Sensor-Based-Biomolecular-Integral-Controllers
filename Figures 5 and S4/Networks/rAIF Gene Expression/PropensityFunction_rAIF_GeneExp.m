function Prop = PropensityFunction_rAIF_GeneExp(x, Parameters)
% Propensity Function for Gene Expression Process Controlled by
% Reference-Based AIF Controller
% 	 Species: 		 X = [X_1; X_2; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	X_1 + X_2		[k_1*X_1]
% 				    R2:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R3:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R4:		 0                  --> 	X_1				[k*Z_1]
% 				    R5:		 0                  --> 	Z_1				[mu]
% 				    R6:		 X_2                --> 	X_2 + Z_2		[theta*X_2]
% 				    R7:		 Z_1 + Z_2          --> 	0				[eta*Z_1*Z_2]
% 				    R8:		 Z_1                --> 	0				[delta*Z_1]
% 				    R9:		 Z_2                --> 	0				[delta*Z_2]

%% Extract Parameters
k_1 = Parameters.k_1;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k = Parameters.k;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
delta = Parameters.delta;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
Z_1 = x(3);
Z_2 = x(4);

%% Propensities
Prop = [ ...
        k_1*X_1; ...
        gamma_1*X_1; ...
		gamma_2*X_2; ...
        k*Z_1; ...
        mu; ...
        theta*X_2; ...
        eta*Z_1*Z_2; ...
        delta*Z_1; ...
        delta*Z_2; ...
	];
end