%%%

% Constants
P = 101325;   % Atmospheric Pressure [Pa]
T = 500;      % Temperature [K]
W1 = 2.0159;  % H2 [g/mol]
W2 = 31.9988; % 02 [g/mol]
W3 = 28.0152; % N2 [g/mol]
R0 = 8314;    %    [kJ/mol*K]

% Initialise Spacial & Temporal Mesh
dx = 1e-6;             
x = 0:dx:1e-4;  % [m]

% Initialise Variable Arrays 
% Yi = Mass Fraction, Xi = molar fraction, W = molar mass of mixture, rhoi = density of species
% Cp = Specific Heat, lamda = Thermal Conductivity, Di = Species Diffusion

Y1 = zeros(length(x),1);
Y2 = zeros(length(x),1);
Y3 = zeros(length(x),1);

X1 = zeros(length(x),1);
X2 = zeros(length(x),1);
X3 = zeros(length(x),1);

W = zeros(length(x),1);
rho_m = zeros(length(x),1);

rho1 = zeros(length(x),1);
rho2 = zeros(length(x),1);
rho3 = zeros(length(x),1);

Cp = zeros(length(x),1);
lambda = zeros(length(x),1);

D1 = zeros(length(x),1);
D2 = zeros(length(x),1);
D3 = zeros(length(x),1);

% Set Initial Solution at t = 0
Y1(1) = 0.4;
Y1(end) = 0;

Y2(1) = 0.4;
Y2(end) = 0;

Y3(1) = 0.2;
Y3(end) = 1;

W(1) = molarmass(Y1(1), Y2(1), Y3(1), W1, W2, W3);
W(end) = molarmass(Y1(end), Y2(end), Y3(end), W1, W2, W3);

X1(1) = W(1) / W1 * Y1(1);
X1(end) = W(end) / W1 * Y1(end);

X2(1) = W(1) / W2 * Y2(1);
X2(end) = W(end) / W2 * Y2(end);

X3(1) = W(1) / W3 * Y3(1);
X3(end) = W(end) / W3 * Y3(end);

X1 = linspace(X1(1), X1(end), length(x));
X2 = linspace(X2(1), X2(end), length(x));
X3 = linspace(X3(1), X3(end), length(x));

figure('Name','Mole Fractions at t=0')
hold on 
plot(x,X1(:), 'DisplayName', 'Hydrogen [H2]')
plot(x,X2(:), 'DisplayName', 'Oxygen [02]')
plot(x,X3(:), 'DisplayName', 'Nitrogen [N2]')
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


