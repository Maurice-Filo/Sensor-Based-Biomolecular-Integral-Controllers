function Prop = PropensityFunction_FP_GeneExpMature(x, Parameters)
% Propensity Function for Gene Expression with Protein Maturation Process Controlled by
% Filtered Proportional Controller
% 	 Species: 		 X = [X_1; X_2; X_3; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	X_1 + X_2		[k_1*X_1]
%                   R2:		 X_2				--> 	X_3		        [c*X_2]
% 				    R3:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R4:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R5:		 X_3				--> 	0				[gamma_3*X_3]
% 				    R6:		 0                  --> 	X_1				[alpha/(1 + Z_2/kappa)]
% 				    R7:		 X_3                --> 	X_3 + Z_2		[theta*X_3]
% 				    R8:		 Z_2                --> 	0				[delta*Z_2]

%% Extract Parameters
k_1 = Parameters.k_1;
c = Parameters.c;
gamma_1 = Parameters.gamma_1;
gamma_2 = Parameters.gamma_2;
gamma_3 = Parameters.gamma_3;
alpha = Parameters.alpha;
kappa = Parameters.kappa;
theta = Parameters.theta;
delta = Parameters.delta;

%% Extract State Variables
X_1 = x(1);
X_2 = x(2);
X_3 = x(3);
Z_2 = x(4);

%% Propensities
Prop = [ ...
        k_1*X_1; ...
        c*X_2; ...
        gamma_1*X_1; ...
		gamma_2*X_2; ...
        gamma_3*X_3; ...
        alpha/(1 + Z_2/kappa); ...
        theta*X_3; ...
        delta*Z_2; ...
	];
end