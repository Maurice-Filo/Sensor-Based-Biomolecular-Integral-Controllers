function Prop = PropensityFunction_P_BirthDeath(x, Parameters)
% Propensity Function for Birth-Death Process Controlled by
% Proportional Controller
% 	 Species: 		 X = [X_1]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		 0                  --> 	X_1				[/(1 + X_1/kappa)]
%% Extract Parameters
gamma_1 = Parameters.gamma_1;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Extract State Variables
X_1 = x(1);

%% Propensities
Prop = [ ...
        gamma_1*X_1; ...
       alpha/(1 + X_1/kappa); ...
	];
end