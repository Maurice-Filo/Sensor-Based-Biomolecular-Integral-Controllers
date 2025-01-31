function Prop = PropensityFunction_P_TwoSpecies(x, Parameters)
% Propensity Function for Two Species Process Controlled by
% Pure Proportional Controller
% 	 Species: 		 X = [X_1; X_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
%                   R2:		 X_1				--> 	X_1 + X_2       [k_1*X_1]
%                   R3:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R4:		 0                  --> 	X_1				[alpha/(1 + X_2/kappa)]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
k_1 = Parameters.k_1;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);

%% Propensities
Prop = [ ...
        gamma_1*X_1; ...
        k_1*X_1; ...
        gamma_2*X_2; ...
        alpha/(1 + X_2/kappa); ...
	];
end