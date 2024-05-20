function lambda_ck = Lambda_CK(Sp,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Sp: Name of the species
%     -- T: Temperature at which the Cp will be evaluated
% OUTPUT:
%     -- Cp_ck: The specific heat capacity at constant pressure for Sp at
%               temperature T.

% Universal gas constant 
R = 8.314510 ; % Unit: J/(mol*K)
Cp = CP_CK(Sp,T);  % Unit: J/(mol*K)             
Cv = Cp - R; % Specific heat under constant volume, % Unit: J/(mol*K)
gama = Cp./Cv; 

Cv = Cv./Sp.M; % Specific heat under constant volume, % Unit: Kg/(mol*K)                 
mu = Mu_CK(Sp,T);

lambda_ck = 0.25.*(9.*gama-5).*mu.*Cv; % Return lambda in W/(m*K) 
