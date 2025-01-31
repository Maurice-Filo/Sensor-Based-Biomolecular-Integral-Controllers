function S = StoichiometryMatrix_FP_GeneExp()
% Stoichiometry Matrix for Gene Expression Process Controlled by
% Filtered Proportional Controller
% 	 Species: 		 X = [X_1; X_2; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	X_1 + X_2		[k_1*X_1]
% 				    R2:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R3:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R4:		 0                  --> 	X_1				[alpha/(1 + Z_2/kappa)]
% 				    R5:		 X_2                --> 	X_2 + Z_2		[theta*X_2]
% 				    R6:		 Z_2                --> 	0				[delta*Z_2]
				    
S = [   0,     -1,      0,      1,      0,      0; ...
        1,      0,     -1,      0,      0,      0; ...
        0,      0,      0,      0,      1,     -1; ...
    ];
end