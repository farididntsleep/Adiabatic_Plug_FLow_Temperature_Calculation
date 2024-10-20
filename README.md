# Reactor Inlet Temperature Calculation with Newton-Raphson and Simpson's Integration

## Description
This MATLAB program calculates the inlet temperature (T0) of an adiabatic plug-flow reactor using the **Newton-Raphson** method for root-finding. It also employs **Simpson's rule** for numerical integration to evaluate the reactor volume.

## How to Use
To run this program, you need to input the required reaction parameters as described below. The program will iteratively solve for the inlet temperature that achieves the desired conversion \( x_{out} \) at the reactor outlet.

### Required Input Parameters
1. **Volumetric flow rate (Fv)** in liters per minute (L/min).
2. **Inlet concentration of A (CA0)** in gram mol per liter (gmol/L).
3. **Density (rho)** in kilograms per liter (kg/L).
4. **Heat capacity (Cp)** in kilocalories per kilogram Kelvin (kcal/kg.K).
5. **Pre-exponential factor (A)** in liters per gram mol per minute (L/gmol/min).
6. **Activation energy (E)** in calories per gram mol (cal/gmol).
7. **Enthalpy change of reaction (delta_Hr)** in kilocalories per gram mol (kcal/gmol).
8. **Reactor volume (Vol)** in liters (L).
9. **Outlet conversion (xout)** between 0 and 1.
10. **Number of integration intervals (n_intervals)** for Simpson's rule.
11. **Initial guess for inlet temperature (T0)** in Kelvin (K).
