function Parameters = Parameters_rAIF_GeneExp()
% Initial Parameters for Gene Expression Process Controlled by
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

Parameters.k_1 = 1;
Parameters.gamma_1 = 1;
Parameters.gamma_2 = 0.1;
Parameters.k = 1;
Parameters.mu = 10;
Parameters.theta = 1;
Parameters.eta = 100;
Parameters.delta = 0;
end

