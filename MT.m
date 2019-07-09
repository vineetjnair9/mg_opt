function [C_MT, Emissions_MT] = MT(P_MT)
% Model for micro gas turbine (MT) - Wu et al. 2014
L = 9.7e-3; % Thermal value of gas (kg/kW)
C_nl = 1; % Price of natural gas ($/m^3)

% Simpler model for fuel consumption cost in controllable DG (MT, DE, FC)
% C_MT = C_nl * (E_GT_total / eta_e);

% Actual source: https://www.epa.gov/sites/production/files/2015-07/documents/catalog_of_chp_technologies_section_5._characterization_-_microturbines.pdf
% Another useful source: https://www.energy.gov/sites/prod/files/2016/09/f33/CHP-Microturbines_0.pdf 

% All values based on HHV of NG
P_rated_MT = 28; % kW - but this is much larger than the peak load
Fuel_MT = 0.434; % Fuel input in MMBtu/h of online operation
C_NG = 2.19; % Cost of natural gas per MMBtu - https://markets.businessinsider.com/commodities/natural-gas-price
eta_e = 0.247; % Electric efficiency
eta_t = 0.469; % Thermal efficiency
eta_overall = 0.716; % Overall efficiency

IC_MT = 3220; % Initial installed capital cost (inc. complete MT package + construction/installation) in $/kW
Fixed_OM_mt = (2/100)*IC_MT;
Variable_OM_MT = 0.013; % $/kWh 

Emissions_MT = 731.6445; % kgCO2(e) / MWh