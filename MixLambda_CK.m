function mixLambda = MixLambda_CK(Mix,X,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Mix: Names of the species in the mixture, e.g. [CO2,H2O,O2,N2]
%     -- X: mole fractions of each species, e.g. [0.2,0.3,0.1,0.4]
%     -- T: Temperature at which the viscosity will be evaluated, in K
% OUTPUT:
%     -- mixLambda: The viscosity for Mix at temperature T.

tLambda_times = 0.0;
tLambda_divide = 0.0;
N = length(Mix);
for a = 1 : N
   
    tLambda_times  = tLambda_times + X(a).*Lambda_CK(Mix(a),T);
    tLambda_divide  = tLambda_divide + X(a)./Lambda_CK(Mix(a),T);
end

mixLambda = 0.5.*(tLambda_times + 1./tLambda_divide); % Return viscosity in Pa.S, Kg/(s m) 


