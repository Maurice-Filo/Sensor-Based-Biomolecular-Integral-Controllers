function Prop = PropensityFunction_rAIF_BirthDeath(x, Parameters)
% Propensity Function for Birth-Death Process Controlled by
% Reference-Based AIF Controller
% 	 Species: 		 X = [X_1; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		 0                  --> 	X_1				[k*Z_1]
% 				    R3:		 0                  --> 	Z_1				[mu]
% 				    R4:		 X_1                --> 	X_1 + Z_2		[theta*X_1]
% 				    R5:		 Z_1 + Z_2          --> 	0				[eta*Z_1*Z_2]
% 				    R6:		 Z_1                --> 	0				[delta*Z_1]
% 				    R7:		 Z_2                --> 	0				[delta*Z_2]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
k = Parameters.k;
mu = Parameters.mu;
theta = Parameters.theta;
eta = Parameters.eta;
delta = Parameters.delta;

%% Extract State Variables
X_1 = x(1);
Z_1 = x(2);
Z_2 = x(3);

%% Propensities
Prop = [ ...
        gamma_1*X_1; ...
        k*Z_1; ...
        mu; ...
        theta*X_1; ...
        eta*Z_1*Z_2; ...
        delta*Z_1; ...
        delta*Z_2; ...
	];
end