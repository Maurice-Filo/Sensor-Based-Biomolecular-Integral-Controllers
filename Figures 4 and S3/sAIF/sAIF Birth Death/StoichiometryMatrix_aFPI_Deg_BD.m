function S = StoichiometryMatrix_aFPI_Deg_BD()
% Stoichiometry Matrix for Network aFPI_Deg_BD
% 	 Species: 		 X = [X_1; Z_1; Z_2]
% 	 Reactions: 	R1:     X_1				--> 	0				[gamma_1*X_1]
% 				    R2:		phi				-->     Z_1				[mu]
% 				    R3:		phi             -->     Z_2				[theta*X_1]
% 				    R4:		Z_1 + Z_2		--> 	0				[eta*Z_1*Z_2]
% 				    R5:		phi				--> 	X_1				[alpha]
% 				    R6:		X_1				--> 	0				[gamma*Z_2 * (X_1 / (kappa_1 + X_1))]

S = [  -1,		0,		0,		0,		1,     -1; ...
		0,		1,		0,	   -1,		0,		0; ...
		0,		0,		1,	   -1,		0,		0];
end