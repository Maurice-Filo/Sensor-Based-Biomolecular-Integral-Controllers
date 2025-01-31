function S = StoichiometryMatrix_P_TwoSpecies()
% Stoichiometry Matrix for Two Species Process Controlled by
% Pure Proportional Controller
% 	 Species: 		 X = [X_1; X_2]
% 	 Reactions: 	R1:		 X_1				--> 	0				[gamma_1*X_1]
%                   R2:		 X_1				--> 	X_1 + X_2       [k_1*X_1]
%                   R3:		 X_2				--> 	0				[gamma_2*X_2]
% 				    R4:		 0                  --> 	X_1				[alpha/(1 + X_2/kappa)]
				    
S = [  -1,      0,      0,      1; ...
        0,      1,     -1,      0; ...
    ];
end