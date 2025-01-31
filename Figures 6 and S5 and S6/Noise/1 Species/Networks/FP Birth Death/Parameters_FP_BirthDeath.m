function Parameters = Parameters_FP_BirthDeath()
% Initial Parameters for Birth-Death Process Controlled by
% Filtered Proportional Controller
% 	 Species: 		 X = [X_1; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		 0                  --> 	X_1				[alpha/(1 + Z_2/kappa)]
% 				    R3:		 X_1                --> 	X_1 + Z_2		[theta*X_1]
% 				    R4:		 Z_2                --> 	0				[delta*Z_2]

Parameters.gamma_1 = 0.1;
Parameters.alpha = 2;
Parameters.kappa = 0.05;
Parameters.theta = 1;
Parameters.delta = 1;
end

