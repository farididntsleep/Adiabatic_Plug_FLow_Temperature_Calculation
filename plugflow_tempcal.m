clc;
clear;
% User input for parameters
Fv = input('Enter the volumetric flow rate (Fv) in L/min: ');
CA0 = input('Enter the inlet concentration of A (CA0) in gmol/L: ');
rho = input('Enter the density (rho) in kg/L: ');
Cp = input('Enter the heat capacity (Cp) in kcal/kg.K: ');
A = input('Enter the pre-exponential factor (A) in L/gmol/min: ');
E = input('Enter the activation energy (E) in cal/gmol: ');
delta_Hr = input('Enter the enthalpy change of the reaction (delta_Hr) in kcal/gmol: ');
R = 1.987;  % Universal gas constant in cal/gmol.K
Vol = input('Enter the reactor volume (Vol) in liters: ');
xout = input('Enter the outlet conversion (xout): ');
n_intervals = input('Enter the number of integration intervals: ');
T0_guess = input('Enter the initial guess for the inlet temperature (T0) in K: ');

% Function to calculate reaction rate constant k
k = @(T) A * exp(-E / (R * T));

% Function f(T0) for use in Newton-Raphson method
f = @(T0) V_simpson(T0, xout, k, n_intervals, CA0, Fv, rho, Cp, delta_Hr) - Vol;

% Newton-Raphson iteration
tol = 1e-6;  % tolerance
max_iter = 100;
T0 = T0_guess;
iter_data = []; % To store iteration data
for iter = 1:max_iter
    fT0 = f(T0);
    dfT0 = (f(T0 + 0.01) - f(T0)) / 0.01; % Numerical derivative approximation
    T0_new = T0 - fT0 / dfT0;
    
    % Store iteration data for display
    iter_data = [iter_data; T0, fT0];
    
    if abs(T0_new - T0) < tol
        break;
    end
    
    T0 = T0_new;
end

% Display iteration results similar to the provided MATLAB output
fprintf('Iteration process:\n');
fprintf('To (K)        f(To)\n');
fprintf('------------------------\n');
for i = 1:size(iter_data, 1)
    fprintf('%-10.4f   %-10.4f\n', iter_data(i, 1), iter_data(i, 2));
end
fprintf('------------------------\n');
fprintf('Calculated reactor inlet temperature: %.4f K\n', T0);

% Function to calculate the reactor volume using Simpson's rule
function V = V_simpson(T0, xout, k, n_intervals, CA0, Fv, rho, Cp, delta_Hr)
    x_values = linspace(0, xout, n_intervals);
    h = (xout - 0) / (n_intervals - 1);
    integral_value = 0;

    for i = 1:n_intervals
        if i == 1 || i == n_intervals
            weight = 1;
        elseif mod(i, 2) == 0
            weight = 4;
        else
            weight = 2;
        end
        
        % Temperature function
        T_current = T0 - (CA0 * delta_Hr) / (rho * Cp) * x_values(i);
        
        % Calculate integral using Simpson's rule
        integral_value = integral_value + weight / (k(T_current) * (1 - x_values(i))^2);
    end
    
    V = (h / 3) * integral_value * (Fv / CA0);
end
