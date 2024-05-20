function mixEnt = MixEnt_CK(Mix,X,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Mix: Names of the species in the mixture, e.g. [CO2,H2O,O2,N2]
%     -- X: mole fractions of each species, e.g. [0.2,0.3,0.1,0.4]
%     -- T: Temperature at which the viscosity will be evaluated, in K
% OUTPUT:
%     -- mixEnt: The enthalpy for Mix at temperature T.

tmixEnt = 0.0;
N = length(Mix);
for a = 1 : N   
    tmixEnt  = tmixEnt + X(a).*Ent_CK(Mix(a),T);  
end

mixEnt = tmixEnt; % Return enthalpy of mixture in J/(mole)


