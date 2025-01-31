function S = StoichiometryMatrix_rAIF_GeneExpMature()
% Stoichiometry Matrix for Gene Expression with Protein Maturation Process Controlled by
% Reference-Based AIF Controller
% 	 Species: 		 X = [X_1; X_2; X_3; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	X_1 + X_2		[k_1*X_1]
%                   R2:		 X_2				--> 	X_3		        [c*X_2]
% 				    R3:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R4:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R5:		 X_3				--> 	0				[gamma_3*X_3]
% 				    R6:		 0                  --> 	X_1				[k*Z_1]
% 				    R7:		 0                  --> 	Z_1				[mu]
% 				    R8:		 X_3                --> 	X_3 + Z_2		[theta*X_3]
% 				    R9:		 Z_1 + Z_2          --> 	0				[eta*Z_1*Z_2]
% 				    R10:	 Z_1                --> 	0				[delta*Z_1]
% 				    R11:	 Z_2                --> 	0				[delta*Z_2]
				    
S = [   0,      0,     -1,      0,      0,      1,      0,      0,      0,      0,      0; ...
        1,     -1,      0,     -1,      0,      0,      0,      0,      0,      0,      0; ...
        0,      1,      0,      0,     -1,      0,      0,      0,      0,      0,      0; ...
        0,      0,      0,      0,      0,      0,      1,      0,     -1,     -1,      0; ...
        0,      0,      0,      0,      0,      0,      0,      1,     -1,      0,     -1; ...
    ];
end