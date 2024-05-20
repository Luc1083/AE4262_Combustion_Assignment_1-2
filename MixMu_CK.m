function mixMu = MixMu_CK(Mix,X,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Mix: Names of the species in the mixture, e.g. [CO2,H2O,O2,N2]
%     -- X: mole fractions of each species, e.g. [0.2,0.3,0.1,0.4]
%     -- T: Temperature at which the viscosity will be evaluated, in K
% OUTPUT:
%     -- mixMu: The viscosity for Mix at temperature T.

tmixMu = 0.0;
N = length(Mix);
for a = 1 : N
    tXPhi = 0.0;
    for b = 1 : N
        tXPhi = tXPhi + X(b).*Phi_AB(Mix(a),Mix(b),T);
    end    
    tmixMu = tmixMu + X(a).*Mu_CK(Mix(a),T)./tXPhi;
end

mixMu = tmixMu; % Return viscosity in Pa.S, Kg/(m*s) 


