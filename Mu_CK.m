function Mu_ck = Mu_CK(Sp,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Sp: Name of the species
%     -- T: Temperature at which the viscosity will be evaluated
% OUTPUT:
%     -- Mu_ck: The viscosity at constant pressure for Sp at
%               temperature T.
Tep = Sp.TpCoes(1);
Sigma = Sp.TpCoes(2);
Omega = 1.147.*(T./Tep).^(-0.145) + (T./Tep + 0.5)^(-2);

tMu_ck = 2.6693e-05.*(Sp.M.*1000.*T)^0.5/(Sigma.^(2).*Omega);

% Cp can be returned in different units, 
% for example, change below to Cp_ck = Cp_J; to return Cp in J/(mole*K)
Mu_ck = 0.1*tMu_ck; % Return Cp in Pa.S, Kg/(s m) 


