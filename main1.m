%%%%%%%%%%%%%%%% Program Information %%%%%%%%%%%%%%%%

% Mixture of 3 species Y1 = H2, Y2 = O2, Y3 = N2
% Posibility of chemical reaction is excluded

% Transport equation used:

% d/dt * rho_i + divergence(rho_i * v) = - divergence(J_i) + M_i * omega_i

% J_i = rho * Y_i * V_i

% Diffusive flux is proportional to gradient of mass fraction

% (rho * Y_i * V_i)_model = - rho * D_i,model * gradient(Y_i)

% determine flux of species having highest mass fraction last

% Sum(Yi * Vi) = 0

% Values to find:
% 1. species mole & mass fractions
% 2. mean molar mass & density
% 3. mixture specific heat & mixture thermal conductivity
% 4. species diffusion coefficients from 4 models

% Initial Conditions. Within domain, mole fractions vary linearly

% x = 0
% Y1 = 0.4
% Y2 = 0.4
% Y3 = 0.2

% x = L = 0.1 mm
% Y1 = 0
% Y2 = 0
% Y3 = 1

%%%%%%%%%%%%%%%% Program Begin %%%%%%%%%%%%%%%%

% Constants
P = 101325;   % Atmospheric Pressure [Pa]
T = 500;      % Temperature [K]
W1 = 2.0159;  % H2 [g/mol]
W2 = 31.9988; % 02 [g/mol]
W3 = 28.0152; % N2 [g/mol]
R0 = 8314;    %    [kJ/mol*K]

% Initialise Spacial & Temporal Mesh
dx = 1e-5;             
dt = 0.1;             
time = 5;               
t = 0:dt:time;  % [s]
x = 0:dx:1e-4;  % [m]

% Initialise Variable Arrays 
% Yi = Mass Fraction, W = molar mass of mixture, Xi = molar fraction, rhoi = density of species i 
Y1 = zeros(length(x),length(t));
Y2 = zeros(length(x),length(t));
Y3 = zeros(length(x),length(t));

W = zeros(length(x),length(t));

X1 = zeros(length(x),length(t));
X2 = zeros(length(x),length(t));
X3 = zeros(length(x),length(t));

rho_m = zeros(length(x),length(t));
rho1 = zeros(length(x),length(t));
rho2 = zeros(length(x),length(t));
rho3 = zeros(length(x),length(t));

% Set Initial Solution at t = 0
Y1(:,1) = linspace(0.4, 0, length(x));
Y2(:,1) = linspace(0.4, 0, length(x));
Y3(:,1) = linspace(0.2, 1, length(x));

W(:,1) = molarmass(Y1(:,1), Y2(:,1), Y3(:,1), W1, W2, W3);
rho_m(:,1) = density(P, R0, T, W(:,1));
rho1(:,1) = rho_m(:,1) .* Y1(:,1);
rho2(:,1) = rho_m(:,1) .* Y2(:,1);
rho3(:,1) = rho_m(:,1) .* Y3(:,1);

X1(:,1) = molarfraction(W1, W(:,1), Y1(:,1));
X2(:,1) = molarfraction(W2, W(:,1), Y2(:,1));
X3(:,1) = molarfraction(W3, W(:,1), Y3(:,1));

% Begin Simulation

% 1. Fick's Diffusion Model, Binary Diffusion Constants (D_i_carrier)

% Hydrogen is Carrier
D_H2_H2 = 33.83 * 10^5;
D_O2_H2 = 18.88 * 10^5;
D_N2_H2 = 17.76 * 10^5;

% Oxygen is Carrier
D_H2_O2 = 18.88 * 10^5;
D_O2_O2 = 5.05 * 10^5;
D_N2_O2 = 4.96 * 10^5;
 
% Nitrogen is Carrier
D_H2_N2 = 17.76 * 10^5;
D_O2_N2 = 4.96 * 10^5;
D_N2_N2 = 4.87 * 10^5;

% 2. Wilke Model

D1_model2 = ((X2 / (D12 * (1-X1))) + (X3 / (D13 * (1-X1)))).^-1;
D2_model2 = ((X1 / (D21 * (1-X2))) + (X3 / (D23 * (1-X2)))).^-1;
D3_model2 = ((X1 / (D31 * (1-X3))) + (X2 / (D32 * (1-X3)))).^-1;

% 3. Le = 1 species diffusivity model (species have same diffusivity as heat)

lambda = MixLambda_CK([H2, O2, N2], [X1(:,i), X2(:,i), X3(:,i)], T);
Cp = MixCp_CK([H2, O2, N2], [X1(:,i), X2(:,i), X3(:,i)], T);
Di_model3 = lambda / (rho * Cp);

% 4. Le = const number model

Le = [0.3, 1.11, 1.0];
lambda = MixLambda_CK([H2, O2, N2], [X1(:,i), X2(:,i), X3(:,i)], T);
Cp = MixCp_CK([H2, O2, N2], [X1(:,i), X2(:,i), X3(:,i)], T);
Di_model4 = lamda / (Le[i] * rho * Cp);

figure('Name','Mole Fractions at t=0')
hold on 
plot(x,X1(:,1), 'DisplayName', 'Hydrogen [H2]')
plot(x,X2(:,1), 'DisplayName', 'Oxygen [02]')
plot(x,X3(:,1), 'DisplayName', 'Nitrogen [N2]')
legend('show')

figure('Name','Mass Fractions at t=0')
hold on 
plot(x,Y1(:,1), 'DisplayName', 'Hydrogen [H2]')
plot(x,Y2(:,1), 'DisplayName', 'Oxygen [02]')
plot(x,Y3(:,1), 'DisplayName', 'Nitrogen [N2]')
legend('show')

figure('Name','Density of species at t=0')
hold on
plot(x,rho1(:,1), 'DisplayName', 'Hydrogen [H2]')
plot(x,rho2(:,1), 'DisplayName', 'Oxygen [02]')
plot(x,rho3(:,1), 'DisplayName', 'Nitrogen [N2]')
legend('show')

function W = molarmass(Y1, Y2, Y3, W1, W2, W3)
    W = (Y1/W1 + Y2/W2 + Y3/W3).^-1;
end

function rho = density(P, R0, T, W)
    rho = P * W / (R0 * T);
end

function X = molarfraction(Wi, W, Yi)
    X = W / Wi .* Yi;
end 