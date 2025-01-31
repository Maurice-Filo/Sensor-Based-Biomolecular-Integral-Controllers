function S = StoichiometryMatrix_P_BirthDeath()
% Stoichiometry Matrix for Birth-Death Process Controlled by
% Proportional Controller
% 	 Species: 		 X = [X_1]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		 0                  --> 	X_1				[alpha/(1 + X_1/kappa)]
				    
S = [  -1,      1; ...
    ];
end