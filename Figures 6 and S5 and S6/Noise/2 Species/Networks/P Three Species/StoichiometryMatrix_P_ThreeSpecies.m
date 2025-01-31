function S = StoichiometryMatrix_P_ThreeSpecies()
% Stoichiometry Matrix for Three Species Process Controlled by
% Pure Proportional Controller
% 	 Species: 		 X = [X_1; X_2; X_3]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
%                   R2:		 X_1				--> 	X_1 + X_2       [k_1*X_1]
%                   R3:		 X_2				--> 	0				[gamma_2*X_2]
%                   R4:		 X_2				--> 	X_3				[c*X_2]
%                   R5:		 X_3				--> 	0				[gamma_3*X_3]
% 				    R6:		 0                  --> 	X_1				[alpha/(1 + X_3/kappa)]
				    
S = [  -1,      0,      0,      0,      0,      1; ...
        0,      1,     -1,     -1,      0,      0; ...
        0,      0,      0,      1,     -1,      0; ...
    ];
end