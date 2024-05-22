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

% D_Model = 1 --- Fick's Diffusion Model Using Carrier Species
% D_Model = 2 --- Wilke Model
% D_Model = 3 --- Le = 1 Species Diffusivity Model (species have same diffusivity as heat)
% D_Model = 4 --- Le = Const Numbers Model (species have constant but different Lewis number)

%%%%%%%%%%%%%%%% Program Begin %%%%%%%%%%%%%%%%

% Program Settings (Select Diffusion Model)
D_Model = 1;
PLOT_LIVE = false; % WARNING: DONT USE WITH LARGE n, NEED TO CLOSE 1 BY 1
PLOT_INITIAL_FINAL = true;

% Constants
P = 101325;   % Atmospheric Pressure [Pa]
T = 500;      % Temperature [K]
W1 = 2.0159;  % H2 [g/mol]
W2 = 31.9988; % 02 [g/mol]
W3 = 28.0152; % N2 [g/mol]
R0 = 8314;    %    [kJ/mol*K]

% Initialise Spacial & Temporal Mesh
dx = 1e-6;             
dt = 1e-19;             
time = 300000*dt;               
t = 0:dt:time;  % [s]
x = 0:dx:1e-4;  % [m]

% Initialise Variable Arrays 
% Yi = Mass Fraction, Xi = molar fraction, W = molar mass of mixture, rhoi = density of species
% Cp = Specific Heat, lamda = Thermal Conductivity, Di = Species Diffusion

Y1 = zeros(length(x),length(t));
Y2 = zeros(length(x),length(t));
Y3 = zeros(length(x),length(t));

X1 = zeros(length(x),length(t));
X2 = zeros(length(x),length(t));
X3 = zeros(length(x),length(t));

W = zeros(length(x),length(t));
rho_m = zeros(length(x),length(t));

rho1 = zeros(length(x),length(t));
rho2 = zeros(length(x),length(t));
rho3 = zeros(length(x),length(t));

Cp = zeros(length(x),length(t));
lambda = zeros(length(x),length(t));

D1 = zeros(length(x),length(t));
D2 = zeros(length(x),length(t));
D3 = zeros(length(x),length(t));

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
n = 1;

% Begin Simulation
while n <= length(t)
    %% State at n
    rho_m(:, n) = rho1(:, n) + rho2(:, n) + rho3(:, n);
    Y1(:, n) = rho1(:, n) ./ rho_m(:, n);
    Y2(:, n) = rho2(:, n) ./ rho_m(:, n);
    Y3(:, n) = rho3(:, n) ./ rho_m(:, n);

    W(:,1) = molarmass(Y1(:,1), Y2(:,1), Y3(:,1), W1, W2, W3);
    X1(:,1) = molarfraction(W1, W(:,1), Y1(:,1));
    X2(:,1) = molarfraction(W2, W(:,1), Y2(:,1));
    X3(:,1) = molarfraction(W3, W(:,1), Y3(:,1));
    %% Live plotting of State
    % Plot live if PLOT_LIVE is enabled
    if PLOT_LIVE
        % Plotting Y1
        subplot(3,1,1);
        plot(x, Y1(:, n), 'r');
        title('Y1 vs. x');
        xlabel('x');
        ylabel('Y1');
        ylim([0 1]);
        
        % Plotting Y2
        subplot(3,1,2);
        plot(x, Y2(:, n), 'g');
        title('Y2 vs. x');
        xlabel('x');
        ylabel('Y2');
        ylim([0 1]);
        
        % Plotting Y3
        subplot(3,1,3);
        plot(x, Y3(:, n), 'b');
        title('Y3 vs. x');
        xlabel('x');
        ylabel('Y3');
        ylim([0 1]);
        
        % Wait for the user to close the plot window
        waitfor(gcf);
        
        % Close the current figure window
        close(gcf);
    end
    %% Diffusion Model
    % 1. Fick's Diffusion Model, Binary Diffusion Constants (D_i_carrier)
    if D_Model == 1
        Carrier = 1;
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
    
            if Carrier == 1
                D1(:, n) = D_H2_H2;
                D2(:, n) = D_O2_H2;
                D3(:, n) = D_N2_H2;
            elseif Carrier == 2
                D1(:, n) = D_H2_O2;
                D2(:, n) = D_O2_O2;
                D3(:, n) = D_N2_O2;
            elseif Carrier == 3
                D1(:, n) = D_H2_N2;
                D2(:, n) = D_O2_N2;
                D3(:, n) = D_N2_N2;
            end
    % 2. Wilke Model, Di_Model2
    elseif D_Model == 2
        D1(:,n) = ((X2(:,n) / (D12 * (1-X1(:,n)))) + (X3(:,n) / (D13 * (1-X1(:,n))))).^-1;
        D2(:,n) = ((X1(:,n) / (D21 * (1-X2(:,n)))) + (X3(:,n) / (D23 * (1-X2(:,n))))).^-1;
        D3(:,n) = ((X1(:,n) / (D31 * (1-X3(:,n)))) + (X2(:,n) / (D32 * (1-X3(:,n))))).^-1;
    
    % 3. Le = 1 species diffusivity model (species have same diffusivity as heat)
    
    % Currently n is the time step counter
    % Currently i us the species counter
    else
        if D_Model == 3
            Le = [1., 1., 1.];
        elseif D_Model == 4
            Le = [0.3, 1.11, 1.0];

        lambda(:,n) = MixLambda_CK({'H2', 'O2', 'N2'}, [X1(:,n), X2(:,n), X3(:,n)], T);
        Cp(:,n) = MixCp_CK({'H2', 'O2', 'N2'}, [X1(:,n), X2(:,n), X3(:,n)], T);
        Di_model3 = lambda(:,n) ./ (rho_m(:,n) .* Cp(:,n)); % element-wise division

        D1(:,n) = Di_model3 * Le(1);
        D2(:,n) = Di_model3 * Le(2);
        D3(:,n) = Di_model3 * Le(3);
        end
    end
    
    %% Update State with Derivatives for next step and Increment n
    Yi_list = [Y1(:, n), Y2(:, n), Y3(:, n)];
    Di_list = [D1(:, n), D2(:, n), D3(:, n)];
    
    for i = [1, 2, 3]
        Yi = Yi_list(:, i);
        Di = Di_list(:, i);
        drhoidt = frstDrv(rho_m(:, n), x) .* Di .* frstDrv(Yi, x) ...
        + rho_m(:, n) .* frstDrv(Di, x) .* frstDrv(Yi, x) ...
        + rho_m(:, n) .* Di .* scndDrv(Yi, x);
        
        if i == 1
            rho1(:, n+1) = rho1(:, n) + dt*drhoidt;
        elseif i == 2
            rho2(:, n+1) = rho2(:, n) + dt*drhoidt;
        elseif i == 3
            rho3(:, n+1) = rho3(:, n) + dt*drhoidt;
        end
    end
    n = n + 1;
end

% Plotting initial and final states if PLOT_INITIAL_FINAL is enabled
if PLOT_INITIAL_FINAL
    % Set the title for the whole window
    figure;
    % Plot initial and final states for Y1
    subplot(3,1,1);
    hold on;
    plot(x, Y1(:, 1), 'r--', 'LineWidth', 1.5); % Initial state for Y1
    plot(x, Y1(:, end), 'r', 'LineWidth', 1.5); % Final state for Y1
    hold off;
    title(['Diffusion Model: ', num2str(D_Model), ', Final Time: ', num2str(time), ', dt: ', num2str(dt), ', dx: ', num2str(dx)]);
    xlabel('x');
    ylabel('Y1');
    legend('Initial', 'Final', 'Location', 'best');
    ylim([0 1]);
    xlim([0 1e-4]);
    
    % Plot initial and final states for Y2
    subplot(3,1,2);
    hold on;
    plot(x, Y2(:, 1), 'g--', 'LineWidth', 1.5); % Initial state for Y2
    plot(x, Y2(:, end), 'g', 'LineWidth', 1.5); % Final state for Y2
    hold off;
    xlabel('x');
    ylabel('Y2');
    legend('Initial', 'Final', 'Location', 'best');
    ylim([0 1]);
    xlim([0 1e-4]);
    
    % Plot initial and final states for Y3
    subplot(3,1,3);
    hold on;
    plot(x, Y3(:, 1), 'b--', 'LineWidth', 1.5); % Initial state for Y3
    plot(x, Y3(:, end), 'b', 'LineWidth', 1.5); % Final state for Y3
    hold off;
    xlabel('x');
    ylabel('Y3');
    legend('Initial', 'Final', 'Location', 'best');
    ylim([0 1]);
    xlim([0 1e-4]);
end


function W = molarmass(Y1, Y2, Y3, W1, W2, W3)
    W = (Y1/W1 + Y2/W2 + Y3/W3).^-1;
end

function rho = density(P, R0, T, W)
    rho = P * W / (R0 * T);
end

function X = molarfraction(Wi, W, Yi)
    X = W / Wi .* Yi;
end 

function d2ydx2 = scndDrv(y, x)
    % secondDerivative calculates the second derivative of y with respect to x
    %
    % Inputs:
    %   x - a vector of x-values
    %   y - a vector of y-values corresponding to the x-values
    %
    % Output:
    %   d2ydx2 - a vector of the second derivative of y with respect to x

    % Check that x and y are vectors of the same length
    if length(x) ~= length(y)
        error('x and y must be vectors of the same length');
    end

    % Number of data points
    n = length(x);
    
    % Initialize the output vector
    d2ydx2 = zeros(size(y));
    
    % Calculate the second derivative using central differences for the interior points
    for i = 2:n-1
        d2ydx2(i) = (y(i+1) - 2*y(i) + y(i-1)) / ((x(i+1) - x(i)) * (x(i) - x(i-1)));
    end
    
    % Use forward difference for the first point
    d2ydx2(1) = (y(3) - 2*y(2) + y(1)) / ((x(2) - x(1))^2);
    
    % Use backward difference for the last point
    d2ydx2(n) = (y(n) - 2*y(n-1) + y(n-2)) / ((x(n) - x(n-1))^2);
end

function dydx = frstDrv(y, x)
    % firstDerivative calculates the first derivative of y with respect to x
    %
    % Inputs:
    %   x - a vector of x-values
    %   y - a vector of y-values corresponding to the x-values
    %
    % Output:
    %   dydx - a vector of the first derivative of y with respect to x

    % Check that x and y are vectors of the same length
    if length(x) ~= length(y)
        error('x and y must be vectors of the same length');
    end

    % Number of data points
    n = length(x);
    
    % Initialize the output vector
    dydx = zeros(size(y));
    
    % Calculate the first derivative using central differences for the interior points
    for i = 2:n-1
        dydx(i) = (y(i+1) - y(i-1)) / (x(i+1) - x(i-1));
    end
    
    % Use forward difference for the first point
    dydx(1) = (y(2) - y(1)) / (x(2) - x(1));
    
    % Use backward difference for the last point
    dydx(n) = (y(n) - y(n-1)) / (x(n) - x(n-1));
end