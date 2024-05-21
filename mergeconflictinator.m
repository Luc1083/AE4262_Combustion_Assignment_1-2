function dydt = timeChangeSpecies(x, y, diff_mode)
    dydx = firstDerivative(x, y);
    d2ydx2 = secondDerivative(x, y);

    if diff_mode == 'Ficks'
        Di 


function d2ydx2 = secondDerivative(x, y)
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

function dydx = firstDerivative(x, y)
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