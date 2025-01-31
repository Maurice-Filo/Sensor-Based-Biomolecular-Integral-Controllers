function Prop = PropensityFunction_P_ThreeSpecies(x, Parameters)
% Propensity Function for Three Species Process Controlled by
% Pure Proportional Controller
% 	 Species: 		 X = [X_1; X_2; X_3]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
%                   R2:		 X_1				--> 	X_1 + X_2       [k_1*X_1]
%                   R3:		 X_2				--> 	0				[gamma_2*X_2]
%                   R4:		 X_2				--> 	X_3				[c*X_2]
%                   R5:		 X_3				--> 	0				[gamma_3*X_3]
% 				    R6:		 0                  --> 	X_1				[alpha/(1 + X_3/kappa)]

%% Extract Parameters
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
gamma_3 = Parameters.gamma_3;
c = Parameters.c;
k_1 = Parameters.k_1;
alpha = Parameters.alpha;
kappa = Parameters.kappa;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
X_3 = x(3);

%% Propensities
Prop = [ ...
        gamma_1*X_1; ...
        k_1*X_1; ...
        gamma_2*X_2; ...
        c*X_2; ...
        gamma_3*X_3; ...
        alpha/(1 + X_3/kappa); ...
	];
end