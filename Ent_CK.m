function H_ck = Ent_CK(Sp,T)
% A small function to evaluate thermodynamic data from chemkin database
% INPUT:
%     -- Sp: Name of the species
%     -- T: Temperature at which the enthalpy will be evaluated (K)
% OUTPUT:
%     -- H_ck: The specific enthalpy (options for different units below)

% Universal gas constant 
R = 8.314510e7 ; % Unit: Chemkin unit, ergs/(mole*K), = 8.314510 J/(mol*K)
                 % 1 erg = 1e-10 kilojoule = 1e-7 joule

 tH = 0.0;
 if (T > Sp.T_range(2))
     for i = 1:5
         tH = tH + Sp.Coes(i).*T.^i./i;  % For upper temperature interval
     end
     tH = tH + Sp.Coes(6);
 else
     for i = 1:5
         tH = tH + Sp.Coes(i+7).*T.^i./i; % Fot lower temperature interval
     end
     tH = tH + Sp.Coes(13);
 end
 
H_ergs = tH .* R; % Unit: ergs/(mole)

H_J = H_ergs .* 1e-07; % Unit: J/(mole)

H_kJkg = H_J ./(1000.*Sp.M); % Unit: kJ/(kg)

H_cal = H_ergs .* 1e-010*238.846; % Unit: cal/(mole), 1 kj/mol = 238.846 cal/mol;

% H can be returned in different units, 
% for example, change below to H_ck = H_J; to return H in J/(mole)
H_ck = H_J; % Return H in J/(mole) 
