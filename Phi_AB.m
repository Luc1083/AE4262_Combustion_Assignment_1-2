function phi_ab = Phi_AB(Sp_a,Sp_b,T)
% The equatiion (6) in "A database of selected transport coefficients for combustion Studies"
% Used for calculation of the dymaic viscosity of a mixture
% INPUT:
%     -- Sp_a: Names of species a
%     -- Sp_a: Names of species b
%     -- T: Temperature at which the viscosity will be evaluated
% OUTPUT:
%     -- Phi_ab 
phi_ab = 8.^(-0.5).*(1+Sp_a.M/Sp_b.M).^(-0.5).*(1+(Mu_CK(Sp_a,T)/Mu_CK(Sp_b,T)).^(0.5).*(1+Sp_b.M/Sp_a.M).^(0.25)).^2;



