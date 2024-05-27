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

Cp = zeros(length(x),1);
lambda = zeros(length(x),1);

% Store Diffusion Model Coefficients
D1_model1 = zeros(length(x),1);
D2_model1 = zeros(length(x),1);
D3_model1 = zeros(length(x),1);

D1_model2 = zeros(length(x),1);
D2_model2 = zeros(length(x),1);
D3_model2 = zeros(length(x),1);

D1_model3 = zeros(length(x),1);
D2_model3 = zeros(length(x),1);
D3_model3 = zeros(length(x),1);

D1_model4 = zeros(length(x),1);
D2_model4 = zeros(length(x),1);
D3_model4 = zeros(length(x),1);




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

% Compute Mass Fraction Along the Domain
A = [X1(1)-1-X1(1)*W1/W3, X1(1)*W1/W2 - X1(1)*W1/W3;
    X2(1)*W2/W1-X2(1)*W2/W3, X2(1)-1-X2(1)*W2/W3];
b = [-X1(1)*W1/W3;
    -X2(1)*W2/W3];
x_result = A \ b;

for i = 1:length(x)
    % Define the coefficient matrix A
    A = [X1(i)-1-X1(i)*W1/W3, X1(i)*W1/W2 - X1(i)*W1/W3;
        X2(i)*W2/W1-X2(i)*W2/W3, X2(i)-1-X2(i)*W2/W3];
    b = [-X1(i)*W1/W3;
        -X2(i)*W2/W3];
    x_result = A \ b;
    
    Y1(i) = x_result(1);
    Y2(i) = x_result(1);
    Y3(i) = 1-x_result(1)-x_result(2);
end

figure('Name','Mass Fractions at t=0')
hold on 
plot(x,Y1(:), 'DisplayName', 'Hydrogen [H2]')
plot(x,Y2(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,Y3(:), 'DisplayName', 'Nitrogen [N2]')
legend('show')

% Compute Mixture Molar Mass along the Domain

W = molarmass(Y1, Y2, Y3, W1, W2, W3);
rho_m = density(P, R0, T, W);
rho1 = rho_m .* Y1;
rho2 = rho_m .* Y2;
rho3 = rho_m .* Y3;

figure('Name','Mixture Molar Mass at t=0')
hold on 
plot(x,W(:), 'DisplayName', 'Mixture Molar Mass')
legend('show')

figure('Name','Mixture Molar Mass & Density at t=0')
hold on 
plot(x,rho_m(:), 'DisplayName', 'Mixture Density')
plot(x,rho1, 'DisplayName', 'Hydrogen [H2]')
plot(x,rho2, 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,rho3, 'DisplayName', 'Nitrogen [N2]')
legend('show')

% 1. Fick's Diffusion Model, Binary Diffusion Constants (D_i_carrier)

dY1_dx = derivative_mass_fraction(Y1, x);
dY2_dx = derivative_mass_fraction(Y2, x);
dY3_dx = derivative_mass_fraction(Y3, x);

% % Hydrogen is Carrier
D_H2_H2 = 33.83 * 10^5;
D_O2_H2 = 18.88 * 10^5;
D_N2_H2 = 17.76 * 10^5;

% % Oxygen is Carrier
D_H2_O2 = 18.88 * 10^5;
D_O2_O2 = 5.05 * 10^5;
D_N2_O2 = 4.96 * 10^5;

% % Nitrogen is Carrier
D_H2_N2 = 17.76 * 10^5;
D_O2_N2 = 4.96 * 10^5;
D_N2_N2 = 4.87 * 10^5;

D1_model1
D2_model1
D3_model1

% 2. Wilke Model, Di_Model2

D1_model2 = (X2 ./ (D_H2_O2 * (1-X1)) + X3 / (D_H2_N2 * (1-X1))).^-1;
D2_model2 = (X1 ./ (D_O2_H2* (1-X2)) + X3 / (D_O2_N2 * (1-X2))).^-1;
D3_model2 = (X1 ./ (D_N2_H2 * (1-X3)) + X2 / (D_N2_O2 * (1-X3))).^-1;

mf_1_model2 = -rho1 .* D1_model2 .* dY1_dx;
mf_2_model2 = -rho2 .* D2_model2 .* dY2_dx;
mf_3_model2 = -rho3 .* D3_model2 .* dY3_dx;

% 3. Le = 1 species diffusivity model (species have same diffusivity as heat)

lambda = MixLambda_CK([H2, O2, N2], [X1, X2, X3], T);
Cp = MixCp_CK([H2, O2, N2], [X1, X2, X3], T);
D1_model3 = lambda / (rho_m * Cp);
D2_model3 = lambda / (rho_m * Cp);
D3_model3 = lambda / (rho_m * Cp);

mf_1_model3 = -rho1 .* D1_model3 .* dY1_dx;
mf_2_model3 = -rho2 .* D2_model3 .* dY2_dx;
mf_3_model3 = -rho3 .* D3_model3 .* dY3_dx;

% 4. Le = const number model

Le = [0.3, 1.11, 1.0];
D1_model4 = lamda / (Le(1) * rho_m * Cp);
D2_model4 = lamda / (Le(2) * rho_m * Cp);
D3_model4 = lamda / (Le(3) * rho_m * Cp);

mf_1_model4 = -rho1 .* D1_model4 .* dY1_dx;
mf_1_model4 = -rho2 .* D2_model4 .* dY2_dx;
mf_3_model4 = -rho3 .* D3_model4 .* dY3_dx;


function W = molarmass(Y1, Y2, Y3, W1, W2, W3)
    W = (Y1/W1 + Y2/W2 + Y3/W3).^-1;
end

function rho = density(P, R0, T, W)
    rho = P * W / (R0 * T);
end

function X = molarfraction(Wi, W, Yi)
    X = W / Wi .* Yi;
end

function deltaY = derivative_mass_fraction(Y,x)
    deltaY = zeros(length(x));
    deltaY(1) = (Y(2) - Y(1)) / (dx);
    for i = 2:length(x)-1
        deltaY(i) = (Y(i+1)-Y(i-1)) / (2 * dx);
    end
    deltaY(end) = (Y(end)-Y(end-1)) / (dx);
end