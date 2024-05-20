clear all
%%%%% BASIC DATA (not to be modified)
%% Temperature range and Polynomial coefficients from Chemkin therm.dat file
% *.Coes[1:7]: a1 to a7 for upper temperature interval (*.T_range(2) to *.T_range(3))
% *.Coes[8:14]: a1 to a7 for lower temperature interval (*.T_range(1) to *.T_range(2))
%CH4
CH4.M = 0.01604; % Kg/mol
CH4.T_range = [300,1000,5000];
CH4.Coes = [1.68347900e+00,1.02372400e-02,-3.87512900e-06,6.78558500e-10,-4.50342300e-14,...
             -1.00807900e+04,9.62339500e+00,7.78741500e-01,1.74766800e-02,-2.78340900e-05,...
             3.04970800e-08,-1.22393100e-11,-9.82522900e+03,1.37221900e+01];
CH4.TpCoes = [141.400,3.746]; % [epsilon/kB,sigma]
         
%C2H6
C2H6.M = 0.03; % Kg/mol
C2H6.T_range = [300,1000,4000];
C2H6.Coes = [4.82593800e+00,1.38404300e-02,-4.55725900e-06,6.72496700e-10,-3.59816100e-14,...
             -1.27177900e+04,-5.23950700e+00,1.46253900e+00,1.54946700e-02,5.78050700e-06,...
             -1.25783200e-08,4.58626700e-12,-1.12391800e+04,1.44322900e+01];
C2H6.TpCoes = [252.3,4.302]; % [epsilon/kB,sigma]

%C3H8
C3H8.M = 0.044; % Kg/mol
C3H8.T_range = [300,1000,5000];
C3H8.Coes = [7.52521700e+00,1.88903400e-02,-6.28392400e-06,9.17937300e-10,-4.81241000e-14,...
             -1.64645500e+04,-1.78439000e+01,8.96920800e-01,2.66898600e-02,5.43142500e-06,...
             -2.12600100e-08,9.24333000e-12,-1.39549200e+04,1.93553300e+01];
C3H8.TpCoes = [266.80,4.982]; % [epsilon/kB,sigma]

%C4H10
C4H10.M = 0.058; % Kg/mol
C4H10.T_range = [300,1500,4000];
C4H10.Coes = [1.99878500e+01,1.03728100e-02,-9.61081800e-07,-9.61081800e-07,8.20282800e-14,...
             -2.62557100e+04,-8.83790700e+01,-2.25661800e+00,5.88173200e-02,-4.52578300e-05,...
             2.03711500e-08,-4.07945800e-12,-1.76023300e+04,3.32959500e+01];
C4H10.TpCoes = [350.900,5.206]; % [epsilon/kB,sigma]

%H2
H2.M = 0.002; % Kg/mol
H2.T_range = [300,1000,5000];
H2.Coes = [2.99142300e+00,7.00064400e-04,-5.63382900e-08,-9.23157800e-12,1.58275200e-15,...
             -8.35034000e+02,-1.35511000e+00,3.29812400e+00,8.24944200e-04,-8.14301500e-07,...
             -9.47543400e-11,4.13487200e-13,-1.01252100e+03,-3.29409400e+00];
H2.TpCoes = [38.000,2.920]; % [epsilon/kB,sigma]
         
%N2
N2.M = 0.028; % Kg/mol
N2.T_range = [300,1000,5000];
N2.Coes = [2.92664000e+00,1.48797700e-03,-5.68476100e-07,1.00970400e-10,-6.75335100e-15,...
             -9.22797700e+02,5.98052800e+00,3.29867700e+00,1.40824000e-03,-3.96322200e-06,...
             5.64151500e-09,-2.44485500e-12,-1.02090000e+03,3.95037200e+00];
N2.TpCoes = [141.400,3.746]; % [epsilon/kB,sigma]
         
%H2O
H2O.M = 0.018; % Kg/mol
H2O.T_range = [300,1000,5000];
H2O.Coes = [2.67214600e+00,3.05629300e-03,-8.73026000e-07,1.20099600e-10,-6.39161800e-15,...
            -2.98992100e+04,6.86281700e+00,3.38684200e+00,3.47498200e-03,-6.35469600e-06,...
            6.96858100e-09,-2.50658800e-12,-3.02081100e+04,2.59023300e+00];
H2O.TpCoes = [572.400,2.605]; % [epsilon/kB,sigma]
        
%CO2
CO2.M = 0.044; % Kg/mol
CO2.T_range = [300,1000,5000];
CO2.Coes = [4.45362300e+00,3.14016900e-03,-1.27841100e-06,2.39399700e-10,-1.66903300e-14,...
            -4.89669600e+04,-9.55395900e-01,2.27572500e+00,9.92207200e-03,-1.04091100e-05,...
            6.86668700e-09,-2.11728000e-12,-4.83731400e+04,1.01884900e+01];
CO2.TpCoes = [244.000,3.763]; % [epsilon/kB,sigma]
        
%O2
O2.M = 0.032; % Kg/mol
O2.T_range = [300,1000,5000];
O2.Coes = [3.69757800e+00,6.13519700e-04,-1.25884200e-07,1.77528100e-11,-1.13643500e-15,...
           -1.23393000e+03,3.18916600e+00,3.21293600e+00,1.12748600e-03,-5.75615000e-07,...
           1.31387700e-09,-8.76855400e-13,-1.00524900e+03,6.03473800e+00];
O2.TpCoes = [107.400,3.458]; % [epsilon/kB,sigma]

%AR
AR.M = 0.040; % Kg/mol
AR.T_range = [300,1000,5000];
AR.Coes = [2.50000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,...    
           -7.45375000e+02, 4.36600100e+00, 2.50000000e+00, 0.00000000e+00, 0.00000000e+00,...    
           0.00000000e+00, 0.00000000e+00, -7.45375000e+02, 4.36600100e+00];                   

AR.TpCoes = [93.3,3.542]; % [epsilon/kB,sigma] taken from Turns

%%% END OF TRANSPORT PROPERTIES DATABASE 
% Constants
% Universal gas constant
R_u = 8.314 % J / (mol K)
% Pressure
Pressure = 100000 % Pa
%% BEGIN OF CASE DESCRIPTION (TO BE MODIFIED BY USER)
%%%%% EXAMPLE 1 CASE OF COMBUSTION from reactants to products AS DETERMINED
%%%%% BY EXCEL TOOL
%%%%% It is assumed that the Excel file has been used to determine
%%%%% the composition of reactants and products for a specified fuel 
%%%%% mixture and air factor.
%%%%% Description of case: repeat some info from Excel input here
%%%%% Fuel mixture: similar to Groningen natural gas
%%%%% Excel input (repeated as a comment here)
%%%%% fuel composition (mole fractions): 
%%%%% 0.85 CH4, 0.05 C2H6, 0.015 C3H8, 0.005 C4H10, 0.08 N2
%%%%% Air factor: lambda = 1.1
%%%%% Part 1: Copied from Excel output (used as input for matlab calc)
%%%%% Composition of reactants and products 
%%%%% (EXCEL CANNOT HANDLE AR; ONL CASE WITH ZERO AR CAN BE HANDLED BY THIS
%%%%% EXAMPLE 1
%FUEL.M = 0.01834 % Fuel mean molar mass from Excel v4 A25
%Mix_reactant = [CH4,C2H6,C3H8,C4H10,O2,N2,H2,AR];  % Name of species in the reactant mixture
%X_reactant = [0.06842,0.00238,0.00051,0.00013,0.17827,0.75029,0.0,0.0];  % From Excel v5 K22-Q22 Mole fraction of species in the mixture
%Mix_prod = [CO2,H2O,O2,N2,AR];  % Name of species in the product mixture
%X_prod = [0.07505,0.14641,0.02965,0.74890,0.0];  % From Excel v5 H37-K37 Mole fraction of species in the product mixture
%NMOLEREACTANT_PERMOLEFUEL = 12.32; % field Z22 in Excel v5, only needed for heat of combustion calc
%NMOLEPRODUCT_PERMOLEFUEL = 12.37; % field U37 in Excel v5 only needed for heat of combustion calc

%%%%% EXAMPLE 2 CASE OF NONREACTING MIXTURE FOR WHICH PROPERTIES ARE TO BE
%%%%% CALCULATED : PURE N2
%Mix_reactant = [CH4,C2H6,C3H8,C4H10,O2,N2,H2,AR];  % Name of species in the reactant mixture
%X_reactant = [0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0];  % Mole fraction of species in the mixture
%Mix_prod = [CO2,H2O,O2,N2,AR];  % Name of species in the product mixture
%X_prod = [0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0];  % Mole fraction of species in the product mixture

%%%%% EXAMPLE 3 CASE OF NONREACTING MIXTURE FOR WHICH PROPERTIES ARE TO BE CALCULATED 
%%%%% H2, O2, AR mixture
Mix_reactant = [CH4,C2H6,C3H8,C4H10,O2,N2,H2,AR];  % Name of species in the reactant mixture
X_reactant = [0.0,0.0,0.0,0.0,0.0349,0.0,0.034603,0.930487];  % Mole fraction of species in the mixture
Mix_prod = [CO2,H2O,O2,N2,AR];  % Name of species in the product mixture
X_prod = [0.0,0.0,0.0,0.0,0.0349,0.0,0.034603,0.930487];  % Mole fraction of species in the product mixture

%%%%% End data 

% Mean molar mass of prod gas
Mixture.M= X_prod(1)*CO2.M + X_prod(2)*H2O.M + X_prod(3)*O2.M + X_prod(4)*N2.M+X_prod(5)*AR.M;
% write prod gas mixture molar mass
display (Mixture.M);

% upper and lower bounds of temperature range
T_prod_min = 300; % Lower bound of temperature range [K]
T_prod_max = 1200; % Upper bound of temperature range [K]

% specified initial temperature of reactants [K]
T_reactant = 298;% specified initial temperature of reactants

%%% END OF CASE DESCRIPTION

%%%%% For testing: Composition of product gas consisting of nitrogen only
%%%%% X_prod = [0.,0.0,0.,1.0];  % Mole fraction of species in the product mixture


%%% START OF SAMPLE CALCULATIONS

%%%% Calculation of product gas properties in the specified temp range
% Calculate thermal conductivity for a range of temperature 
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Lambda(j) = MixLambda_CK(Mix_prod,X_prod,T_flue(j));
end
 
% Calculate dynamic viscosity for a range of temperature 
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Mu(j) = MixMu_CK(Mix_prod,X_prod,T_flue(j));
end

   % Calculate specific heat per mole for a range of temperature 
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Cp(j) = MixCp_CK(Mix_prod,X_prod,T_flue(j));
end   

% Calculate specific heat per kg for a range of temperature 
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Specific_Cp(j) = Cp(j)/Mixture.M ;
end
% Calculate Prandtl number for a range of temperature 
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Pr(j) = MixPr_CK(Mix_prod,X_prod,T_flue(j));
end

% Calculate enthalpy for a range of temperature 
% Enthalpy of formation
T_ref = 298;       
Ent_ref = MixEnt_CK(Mix_prod,X_prod,T_ref);
Specific_Ent_ref = Ent_ref/Mixture.M ;
% Enthalpy (sum of enthalpy of formation and sensible enthalpy)
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Ent(j) = MixEnt_CK(Mix_prod,X_prod,T_flue(j));
       Specific_Ent(j) = Ent(j)/Mixture.M ;
end
% Sensible enthalpy
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Ent_Sens(j) = Ent(j)-Ent_ref;
       Specific_Ent_Sens(j) = Ent_Sens(j)/Mixture.M ;
end

% Calculate density for a range of temperature 
T_flue = T_prod_min:100:T_prod_max;
for j =1:length(T_flue)
       Rho(j) = (Mixture.M * Pressure) / (R_u * T_flue(j));
end

% Plot results
figure (1)
 plot(T_flue,Lambda)
 xlabel('Temperature [K]')
 ylabel('Lambda flue [W/(m*K)]')
figure (2)
 plot(T_flue, Mu)
 xlabel('Temperature [K]')
 ylabel('Dynamic viscosity [Pa s]')
figure (3)
 plot(T_flue, Cp)
 xlabel('Temperature [K]')
 ylabel('Specific heat per mole [(J/(mole * K)]')
 figure (4)
 plot(T_flue, Specific_Cp)
 xlabel('Temperature [K]')
 ylabel('Specific heat per kg[(J/(kg * K)]')
figure (5) 
 plot(T_flue, Pr)
 xlabel('Temperature [K]')
 ylabel('Prandt number [-]')
figure (6) 
 plot(T_flue, Ent)
 xlabel('Temperature [K]')
 ylabel('Enthalpy [(J/(mole)]')
 figure (7) 
 plot(T_flue, Specific_Ent)
 xlabel('Temperature [K]')
 ylabel('Specific Enthalpy [(J/(kg)]')
 figure (8) 
 plot(T_flue, Ent_Sens)
 xlabel('Temperature [K]')
 ylabel('Sensible Enthalpy [(J/(mole)]')
 figure (9) 
 plot(T_flue, Specific_Ent_Sens)
 xlabel('Temperature [K]')
 ylabel('Specific Sensible Enthalpy [(J/(kg)]')
 figure (10) 
 plot(T_flue, Rho)
 xlabel('Temperature [K]')
 ylabel('Mass density [(kg/m3)]')
 figure (11)
 plot(T_flue, Specific_Cp)
 xlabel('Temperature [K]')
 ylabel('Specific heat per kg[(J/(kg * K)]')
 
%%  
%% ADIABATIC FLAME TEMPERATURE CALCULATION (NEEDS CORRECT MIXTURE OF REACTANTS AND PRODUCTS (AS IN EXAMPLE1)
%% FOR EXAMPLES 2 AND 3 IT CANNOT BE USED
%% Using MixEnt_CK and an iterative procedure to calculate temperature number of a products
%% Input: composition and temperature of reactants
%% Input: composition of products
%% Output: temperature of products (adiabatic flame temperature)
%T_prod = 1000;               % A guessed initial temperature for the combustion products
%T_uplim = 3000;              % upper limit of products temperature
%T_lowlim = T_reactant;              % low limit of products temperature
%Ent_reactant = MixEnt_CK(Mix_reactant,X_reactant,T_reactant); % Calculate the enthalpy of reactants
%Ent_prod = MixEnt_CK(Mix_prod,X_prod,T_prod); % Calculate the enthalpy of products
%Ent_tol = 1.0; % Tolerance for the difference of enthalpy between reactants and products 
%while (abs(Ent_prod - Ent_reactant) > Ent_tol)
%    if (Ent_prod > Ent_reactant)
%        T_uplim = T_prod;
%        T_prod = 0.5*(T_prod + T_lowlim);
%    else
%        T_lowlim = T_prod;
%        T_prod = 0.5*(T_prod + T_uplim);
%    end
%    Ent_prod = MixEnt_CK(Mix_prod,X_prod,T_prod); % Calculate the enthalpy of products
%end
%T_adiabatic = T_prod
%% The adiabatic temperature of products
%display(T_adiabatic)
