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
% 5. 

% Initial Conditions. Within domain, mole fractions vary linearly

% x = 0
% Y1 = 0.4
% Y2 = 0.4
% Y3 = 0.2

% x = L = 0.1 mm
% Y1 = 0
% Y2 = 0
% Y3 = 1

% Constants
p_0 = 101325; % Atmospheric Pressure [Pa]
T_0 = 500;    % Temperature [K]

% Initial Solution
dx = 0.01;
x = 0:dx:0.1;
Y1 = linspace(0.4, 0, length(x));
Y2 = linspace(0.4, 0, length(x));
Y3 = linspace(0.2, 1, length(x));

figure('Name','Mole Fractions at t=0')
hold on 
plot(x,Y1)
plot(x,Y2)
plot(x,Y3)