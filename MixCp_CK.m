function mixCp = MixCp_CK(Mix,X,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Mix: Names of the species in the mixture, e.g. [CO2,H2O,O2,N2]
%     -- X: mole fractions of each species, e.g. [0.2,0.3,0.1,0.4]
%     -- T: Temperature at which the viscosity will be evaluated, in K
% OUTPUT:
%     -- mixLambda: The specific heat for Mix at temperature T.

tmixCp = 0.0;

N = length(Mix);
mixM = 0.0;
for a = 1 : N   
    tmixCp  = tmixCp + X(a).*CP_CK(Mix(a),T);  
    mixM = mixM + X(a).*Mix(a).M;
end

mixCp = tmixCp./mixM; % Return Cp of mixture in J/(kg*K)
% mixCp = tmixCp; % Return Cp of mixture in J/(mole*K)


