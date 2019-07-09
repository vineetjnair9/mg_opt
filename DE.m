%% Diesel generator sub-model

% Diesel generator specs from Moshi et al. 2016
P_DE_rated = 16; % Rated maximum power of diesel generator (kW) 
P_DE_min = 4.8; % Lower limit on power (kW) or k_gen * P_rated (k_gen = 0.3) from Fathima et al. 
L_DE = 15000; % Lifetime (h)

% Fuel consumption will also depend on overall efficiency of system
eta_brake_thermal = 0.5; % Brake thermal efficiency of diesel engine
eta_generator = 0.9; % Conversion efficiency of electric generator = electric efficiency
eta_overall = eta_brake_thermal * eta_generator; % From Fathima et al. 

% Fuel consumption costs as a quadratic function of power output (Nemati et al.)

% General model for all thermal generators (DE and MT) - Xiao et al. 2014
% F = F0 * P_rated + F1 * P - F0,F1 are coefficients of fuel % consumption-power curve
fuel_DE = 0.06 * P_DE_rated + 0.0246 * P_DE; % fuel consumption (L/h) where P is in kW 
% Diesel costs ($) assuming a price of $3.20/US liquid gallon (value for Tanzania from NREL ReOpt)
C_fuel_DE = 3.20 * (fuel_DE / 3.78541); % Diesel fuel cost per h 

IC_DE = 12500; % Initial capital cost ($)
RC_DE = 11000; % Replacement cost ($)
%IC_DE = 32000; % Initial capital cost ($) - https://dieselgeneratordirect.uk/
Fixed_OM_DE = (2/100)*IC_DE; % Fixed O&M costs ($/y)
Variable_OM_DE = 0.24; % Variable O&M costs ($/h of operation online)
SUC_DE = 0.45; % Start-up cost ($) - not dependent on the amount of power produced?
SDC_DE = 0.23; % Shut-down cost ($)

% Emissions - DE's GHG emissions factors (kg of CO2(e)/kWh)
Emissions_DE = 3.5; % kg of CO2 / L of diesel fuel consumed - R. Dufo-López et al. 2011
