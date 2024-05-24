4function Cp_ck = CP_CK(Sp,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Sp: Name of the species
%     -- T: Temperature at which the Cp will be evaluated
% OUTPUT:
%     -- Cp_ck: The specific heat capacity at constant pressure for Sp at
%               temperature T.

% Universal gas constant 
R = 8.314510e7 ; % Unit: Chemkin unit, ergs/(mole*K), = 8.314510 J/(mol*K)
                 % 1 erg = 1e-10 kilojoule = 1e-7 joule

 tCp = 0.0;
 if (T > Sp.T_range(2))
     for i = 1:5
         tCp = tCp + Sp.Coes(i).*T.^(i-1);  % For upper temperature interval
     end
 else
     for i = 1:5
         tCp = tCp + Sp.Coes(i+7).*T.^(i-1); % Fot lower temperature interval
     end
 end

Cp_ergs = tCp .* R; % Unit: ergs/(mole*K)

Cp_J = Cp_ergs .* 1e-07; % Unit: J/(mole*K)

Cp_kJkg = Cp_J ./(1000.*Sp.M); % Unit: kJ/(kg*K)

Cp_cal = Cp_ergs .* 1e-010*238.846; % Unit: cal/(mole*K), 1 kj/mol = 238.846 cal/mol;

% Cp can be returned in different units, 
% for example, change below to Cp_ck = Cp_J; to return Cp in J/(mole*K)
Cp_ck = Cp_J; % Return Cp in J/(mole K) 
