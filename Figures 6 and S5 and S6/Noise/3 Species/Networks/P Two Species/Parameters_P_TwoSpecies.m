function Parameters = Parameters_P_TwoSpecies()
% Initial Parameters for Two Species Process Controlled by
% Pure Proportional Controller
% 	 Species: 		 X = [X_1; X_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
%                   R2:		 X_1				--> 	X_1 + X_2       [k_1*X_1]
%                   R3:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R4:		 0                  --> 	X_1				[alpha/(1 + X_2/kappa)]

Parameters.gamma_1 = 1;
Parameters.gamma_2 = 0.1;
Parameters.k_1 = 1;
Parameters.alpha = 2;
Parameters.kappa = 0.05;
end

