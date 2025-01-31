function S = StoichiometryMatrix_rAIF_GeneExp()
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
				    
S = [   0,     -1,      0,      1,      0,      0,      0,      0,      0; ...
        1,      0,     -1,      0,      0,      0,      0,      0,      0; ...
        0,      0,      0,      0,      1,      0,     -1,     -1,      0; ...
        0,      0,      0,      0,      0,      1,     -1,      0,     -1; ...
    ];
end