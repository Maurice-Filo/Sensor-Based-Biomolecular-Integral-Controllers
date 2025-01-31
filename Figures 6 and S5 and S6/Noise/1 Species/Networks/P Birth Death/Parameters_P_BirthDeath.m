function Parameters = Parameters_P_BirthDeath()
% Initial Parameters for Birth-Death Process Controlled by
% Proportional Controller
% 	 Species: 		 X = [X_1]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		 0                  --> 	X_1				[alpha/(1 + X_1/kappa)]

Parameters.gamma_1 = 0.1;
Parameters.alpha = 2;
Parameters.kappa = 0.05;
end

