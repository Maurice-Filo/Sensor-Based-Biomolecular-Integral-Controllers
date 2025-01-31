function S = StoichiometryMatrix_rAIF_BirthDeath()
% Stoichiometry Matrix for Birth-Death Process Controlled by
% Reference-Based AIF Controller
% 	 Species: 		 X = [X_1; Z_1; Z_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		 0                  --> 	X_1				[k*Z_1]
% 				    R3:		 0                  --> 	Z_1				[mu]
% 				    R4:		 X_1                --> 	X_1 + Z_2		[theta*X_1]
% 				    R5:		 Z_1 + Z_2          --> 	0				[eta*Z_1*Z_2]
% 				    R6:		 Z_1                --> 	0				[delta*Z_1]
% 				    R7:		 Z_2                --> 	0				[delta*Z_2]
				    
S = [  -1,      1,      0,      0,      0,      0,      0; ...
        0,      0,      1,      0,     -1,     -1,      0; ...
        0,      0,      0,      1,     -1,      0,     -1; ...
    ];
end