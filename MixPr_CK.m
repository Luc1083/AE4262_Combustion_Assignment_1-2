function mixPr = MixPr_CK(Mix,X,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Mix: Names of the species in the mixture, e.g. [CO2,H2O,O2,N2]
%     -- X: mole fractions of each species, e.g. [0.2,0.3,0.1,0.4]
%     -- T: Temperature at which the viscosity will be evaluated, in K
% OUTPUT:
%     -- mixPr: The Prandtl number for Mix at temperature T.

% calculate mixture molar mass
N = length(Mix);
mixM = 0.0;
for a = 1 : N   
    mixM = mixM + X(a).*Mix(a).M;
end

Mu_mix = MixMu_CK(Mix,X,T);      % Call the function
Cp_mix = MixCp_CK(Mix,X,T);      % Specific heat per mole

Lambda_mix = MixLambda_CK(Mix,X,T);
mixPr = Mu_mix * (Cp_mix/mixM)/Lambda_mix; % Return Prandtl number



