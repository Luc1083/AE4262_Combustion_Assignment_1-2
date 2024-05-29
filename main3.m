%% Constants
P = 101325;   % Atmospheric Pressure [Pa]
T = 500;      % Temperature [K]
W1 = 2.0159 / 1000;  % H2 [kg/mol]
W2 = 31.9988 / 1000; % 02 [kg/mol]
W3 = 28.0152 / 1000; % N2 [kg/mol]
W_CH4 = 16.0428 / 1000; % CH4 [kg/mol]
R0 = 8.314;    %    [J/mol*K]

%% Initialise Spacial & Temporal Mesh
dx = 1e-6;             
x = 0:dx:1e-4;  % [m]

%% Initialise Variable Arrays 
% Yi = Mass Fraction, Xi = molar fraction, W = molar mass of mixture

Y1 = zeros(1,length(x));
Y2 = zeros(1,length(x));
Y3 = zeros(1,length(x));

X1 = zeros(1,length(x));
X2 = zeros(1,length(x));
X3 = zeros(1,length(x));

W = zeros(1,length(x));

%% Set Initial Solution at t = 0
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
plot(x,X1(:), 'DisplayName', 'Hydrogen [H2]', 'Color', 'r')
plot(x,X2(:), 'DisplayName', 'Oxygen [02]', 'Color', 'k')
plot(x,X3(:), 'DisplayName', 'Nitrogen [N2]', 'Color', 'b')
xlabel('Domain Position [m]')
ylabel('Mole Fraction [-]')
legend('show')

% Compute Mass Fraction Along the Domain
for i = 1:length(x)
    
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
plot(x,Y1(:), 'DisplayName', 'Hydrogen [H2]', 'Color', 'r')
plot(x,Y2(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,Y3(:), 'DisplayName', 'Nitrogen [N2]', 'Color', 'b')
xlabel('Domain Position [m]')
ylabel('Mass Fraction [-]')
legend('show')

%% Compute Mixture Molar Mass along the Domain

W = molarmass(Y1, Y2, Y3, W1, W2, W3);
rho_m = density(P, R0, T, W);
rho1 = rho_m .* Y1;
rho2 = rho_m .* Y2;
rho3 = rho_m .* Y3;

figure('Name','Mixture Molar Mass at t=0')
hold on 
plot(x,W(:), 'DisplayName', 'Mixture Molar Mass')
xlabel('Domain Position [m]')
ylabel('Molar Mass [kg/mol]')
legend('show')

figure('Name','Mixture & Species Density at t=0')
hold on 
plot(x,rho_m(:), 'DisplayName', 'Mixture Density')
plot(x,rho1, 'DisplayName', 'Hydrogen [H2]', 'Color', 'r')
plot(x,rho2, 'DisplayName', 'Oxygen [02]','LineStyle','--', 'Color','k')
plot(x,rho3, 'DisplayName', 'Nitrogen [N2]','Color','b')
xlabel('Domain Position [m]')
ylabel('Density [kg/m^3]')
legend('show')

%% Compute Molar Concentration along the Domain

M1 = rho_m .* Y1 / W1;
M2 = rho_m .* Y2 / W2;
M3 = rho_m .* Y3 / W3;

figure('Name','Molar Concentration at t=0')
hold on 
plot(x,M1(:), 'DisplayName', 'Hydrogen [H2]','Color', 'r')
plot(x,M2(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,M3(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
xlabel('Domain Position [m]')
ylabel('Molar Concentration [mol/m^3]')
legend('show')

%% Compute Mixture specific heat and mixture thermal conductivity
lambda = zeros(1,length(x));
Cp = zeros(1,length(x));

for i = 1:length(x)
    lambda(i) = MixLambda_CK([H2, O2, N2], [X1(i), X2(i), X3(i)], T);
    Cp(i) = MixCp_CK([H2, O2, N2], [X1(i), X2(i), X3(i)], T);
end

% Create a new figure
figure;

% Plot Cp with primary y-axis
yyaxis left;
plot(x, Cp, 'b', 'LineWidth', 2);
ylabel('Mixture Specific Heat [J/(kg*K)]');
xlabel('X-axis');

% Add secondary y-axis for lambda
yyaxis right;
plot(x, lambda, 'LineWidth', 2);
ylabel('Mixture Conductivity [W/(m*K)]');

% Add legend
legend('Cp', '\lambda', 'Location', 'best');

%% Diffusion Model Values

dY1_dx = derivative_mass_fraction(Y1, x);
dY2_dx = derivative_mass_fraction(Y2, x);
dY3_dx = derivative_mass_fraction(Y3, x);

% % Hydrogen is Carrier
D_H2_H2 = 33.83 / 10^5;
D_O2_H2 = 18.88 / 10^5;
D_N2_H2 = 17.76 / 10^5;
D_CH4_H2 = 17.18 / 10^5;

% % Oxygen is Carrier
D_H2_O2 = 18.88 / 10^5;
D_O2_O2 = 5.05 / 10^5;
D_N2_O2 = 4.96 / 10^5;
D_CH4_O2 = 5.49 / 10^5;

% % Nitrogen is Carrier
D_H2_N2 = 17.76 / 10^5;
D_O2_N2 = 4.96 / 10^5;
D_N2_N2 = 4.87 / 10^5;
D_CH4_N2 = 5.37 / 10^5;

% % Methane is Carrier
D_H2_CH4 = 17.18 / 10^5;
D_O2_CH4 = 5.49 / 10^5;
D_N2_CH4 = 5.37 / 10^5;
D_CH4_CH4 = 5.63 / 10^5;

%% 1. Fick's Diffusion Model, Binary Diffusion Constants (D_i_carrier)

D1_model1 = D_H2_N2 * ones(1,length(x));
D2_model1 = D_O2_N2 * ones(1,length(x));
D3_model1 = D_N2_N2 * ones(1,length(x));

mf_1_model1 = -rho_m .* D1_model1 .* dY1_dx;
mf_2_model1 = -rho_m .* D2_model1 .* dY2_dx;
mf_3_model1 = -mf_2_model1 -mf_1_model1;
% mf_3_model1_nc = -rho_m .* D3_model1 .* dY3_dx;

figure('Name','Diffusion Constants (Fick Model) at t=0')
hold on 
plot(x,D1_model1(:), 'DisplayName', 'Hydrogen [H2]','Color', 'r')
plot(x,D2_model1(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,D3_model1(:), 'DisplayName', 'Nitrogen [N2]','Color', 'b')
xlabel('Domain Position [m]')
ylabel('Diffusion Coefficient [m^2/s]')
legend('show')

figure('Name','Mass Flux (Fick Model) at t=0')
hold on 
plot(x,mf_1_model1(:), 'DisplayName', 'Hydrogen [H2]','Color', 'r')
plot(x,mf_2_model1(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,mf_3_model1(:), 'DisplayName', 'Nitrogen [N2]','Color', 'b')
% plot(x,mf_3_model1_nc(:), 'DisplayName', 'Nitrogen [N2] Non Consistent')
xlabel('Domain Postition [m]')
ylabel('Mass Flux [kg/m^2*s]')
legend('show')

%% 2. Wilke Model, Di_Model2

D1_model2 = (X2 ./ (D_H2_O2 .* (1-X1)) + X3 ./ (D_H2_N2 .* (1-X1))).^-1;
D2_model2 = (X1 ./ (D_O2_H2 .* (1-X2)) + X3 ./ (D_O2_N2 .* (1-X2))).^-1;
D3_model2 = (X1 ./ (D_N2_H2 .* (1-X3)) + X2 ./ (D_N2_O2 .* (1-X3))).^-1;

mf_1_model2 = -rho_m .* D1_model2 .* dY1_dx;
mf_2_model2 = -rho_m .* D2_model2 .* dY2_dx;
mf_3_model2 = -mf_2_model2 -mf_1_model2;
% mf_3_model2_nc = -rho_m .* D3_model2 .* dY3_dx;

figure('Name','Diffusion Constants (Wilke Model) at t=0')
hold on 
plot(x,D1_model2(:), 'DisplayName', 'Hydrogen [H2]','Color','r')
plot(x,D2_model2(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,D3_model2(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
xlabel('Domain Position [m]')
ylabel('Diffusion Coefficient [m^2/s]')
legend('show')

figure('Name','Mass Flux (Wilke Model) at t=0')
hold on 
plot(x,mf_1_model2(:), 'DisplayName', 'Hydrogen [H2]','Color','r')
plot(x,mf_2_model2(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,mf_3_model2(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
% plot(x,mf_3_model2_nc(:), 'DisplayName', 'Nitrogen [N2] Non Consistent')
xlabel('Domain Postition [m]')
ylabel('Mass Flux [kg/m^2*s]')
legend('show')

%% 3. Le = 1 species diffusivity model (species have same diffusivity as heat)

D1_model3 = lambda ./ (rho_m .* Cp);
D2_model3 = lambda ./ (rho_m .* Cp);
D3_model3 = lambda ./ (rho_m .* Cp);

mf_1_model3 = -rho_m .* D1_model3 .* dY1_dx;
mf_2_model3 = -rho_m .* D2_model3 .* dY2_dx;
mf_3_model3 = -mf_2_model3 -mf_1_model3;
% mf_3_model3_nc = -rho_m .* D3_model3 .* dY3_dx;

figure('Name','Diffusion Constants (Le = 1) at t=0')
hold on 
plot(x,D1_model3(:), 'DisplayName', 'Hydrogen [H2]','Color','r')
plot(x,D2_model3(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,D3_model3(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
xlabel('Domain Position [m]')
ylabel('Diffusion Coefficient [m^2/s]')
legend('show')

figure('Name','Mass Flux (Le = 1) at t=0')
hold on 
plot(x,mf_1_model3(:), 'DisplayName', 'Hydrogen [H2]','Color','r')
plot(x,mf_2_model3(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,mf_3_model3(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
% plot(x,mf_3_model3_nc(:), 'DisplayName', 'Nitrogen [N2] Non Consistent', 'LineStyle','--')
xlabel('Domain Postition [m]')
ylabel('Mass Flux [kg/m^2*s]')
legend('show')

%% 4. Le = const number model

D1_model4 = lambda ./ (0.3 .* rho_m .* Cp);
D2_model4 = lambda ./ (1.11 .* rho_m .* Cp);
D3_model4 = lambda ./ (1.0 .* rho_m .* Cp);

mf_1_model4 = -rho_m .* D1_model4 .* dY1_dx;
mf_2_model4 = -rho_m .* D2_model4 .* dY2_dx;
mf_3_model4 = -mf_2_model4 -mf_1_model4;
mf_3_model4_nc = -rho_m .* D3_model4 .* dY3_dx;

figure('Name','Diffusion Constants (Le = const) at t=0')
hold on 
plot(x,D1_model4(:), 'DisplayName', 'Hydrogen [H2]','Color','r')
plot(x,D2_model4(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,D3_model4(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
xlabel('Domain Position [m]')
ylabel('Diffusion Coefficient [m^2/s]')
legend('show')

figure('Name','Mass Flux (Le = const) at t=0')
hold on 
plot(x,mf_1_model4(:), 'DisplayName', 'Hydrogen [H2]','Color','r')
plot(x,mf_2_model4(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,mf_3_model4(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
% plot(x,mf_3_model4_nc(:), 'DisplayName', 'Nitrogen [N2] Non Consistent')
xlabel('Domain Postition [m]')
ylabel('Mass Flux [kg/m^2*s]')
legend('show')

figure('Name','Mass Fraction Gradient dydx')
hold on 
plot(x,dY1_dx(:), 'DisplayName', 'Hydrogen [H2]','Color','r')
plot(x,dY2_dx(:), 'DisplayName', 'Oxygen [02]','LineStyle','--','Color','k')
plot(x,dY3_dx(:), 'DisplayName', 'Nitrogen [N2]','Color','b')
legend('show')

figure('Name','Specific Heat')
hold on 
plot(x,Cp(:))
legend('show')

figure('Name','Thermal Conductivity')
hold on 
plot(x,lambda(:))
legend('show')

%% 3 plots next to each other showing mass flux of diff models
figure('Name','Species Mass Flux')

subplot(1, 3, 1);
hold on
title('Hydrogen (H2)')
plot(x, mf_1_model1(:), 'DisplayName', 'Fick Model')
plot(x, mf_1_model2(:), 'DisplayName', 'Wilke Model')
plot(x, mf_1_model3(:), 'DisplayName', 'Le = 1')
plot(x, mf_1_model4(:), 'DisplayName', 'Le = 0.3')
xlabel('Pos. [m]')
ylabel('Mass Flux [kg/m^2*s]')
legend('show')

subplot(1, 3, 2);
hold on
title('Oxygen (O2)')
plot(x, mf_2_model1(:), 'DisplayName', 'Fick Model')
plot(x, mf_2_model2(:), 'DisplayName', 'Wilke Model')
plot(x, mf_2_model3(:), 'DisplayName', 'Le = 1')
plot(x, mf_2_model4(:), 'DisplayName', 'Le = 1.11')
xlabel('Pos. [m]')
ylabel('Mass Flux [kg/m^2*s]')
legend('show')

subplot(1, 3, 3);
hold on
title('Nitrogen (N2)')
plot(x, mf_3_model1(:), 'DisplayName', 'Fick Model')
plot(x, mf_3_model2(:), 'DisplayName', 'Wilke Model')
plot(x, mf_3_model3(:), 'DisplayName', 'Le = 1')
plot(x, mf_3_model4(:), 'DisplayName', 'Le = 1.0')
xlabel('Pos. [m]')
ylabel('Mass Flux [kg/m^2*s]')
legend('show')

% 3 plots next to each other showing diffusion coefficients of diff models
figure('Name','Species Diffusion Coefficients')

subplot(1, 3, 1);
hold on
title('Hydrogen (H2)')
plot(x, D1_model1(:), 'DisplayName', 'Fick Model')
plot(x, D1_model2(:), 'DisplayName', 'Wilke Model')
plot(x, D1_model3(:), 'DisplayName', 'Le = 1')
plot(x, D1_model4(:), 'DisplayName', 'Le = 0.3')
xlabel('Pos. [m]')
ylabel('Diffusion Coefficient [m^2/s]')
legend('show')

subplot(1, 3, 2);
hold on
title('Oxygen (O2)')
plot(x, D2_model1(:), 'DisplayName', 'Fick Model')
plot(x, D2_model2(:), 'DisplayName', 'Wilke Model')
plot(x, D2_model3(:), 'DisplayName', 'Le = 1')
plot(x, D2_model4(:), 'DisplayName', 'Le = 1.11')
xlabel('Pos. [m]')
ylabel('Diffusion Coefficient [m^2/s]')
legend('show')

subplot(1, 3, 3);
hold on
title('Nitrogen (N2)')
plot(x, D3_model1(:), 'DisplayName', 'Fick Model')
plot(x, D3_model2(:), 'DisplayName', 'Wilke Model')
plot(x, D3_model3(:), 'DisplayName', 'Le = 1')
plot(x, D3_model4(:), 'DisplayName', 'Le = 1.0')
xlabel('Pos. [m]')
ylabel('Diffusion Coefficient [m^2/s]')
legend('show')

function W = molarmass(Y1, Y2, Y3, W1, W2, W3)
    W = (Y1/W1 + Y2/W2 + Y3/W3).^-1;
end

function rho = density(P, R0, T, W)
    rho = P * W / (R0 * T);
end

function deltaY = derivative_mass_fraction(Y, x)
    deltaY = zeros(1,length(x));
    deltaY(1) = (-Y(3) + 4 * Y(2) - 3 * Y(1)) / (2 * (x(2) - x(1)));
    for i = 2:length(x)-1
        deltaY(i) = (Y(i+1)-Y(i-1)) / (x(i+1) - x(i-1));
    end
    deltaY(end) = (3 * Y(end) - 4 * Y(end-1) + Y(end-2)) / (2 * (x(end) - x(end-1)));
end