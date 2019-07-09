   
%% RENEWABLE DISTRIBUTED GENERATORS
% Solar
% Latitude and longitude of microgrid location as user inputs
A_panel = 1.625 * 1.019; % Surface area of 1 panel rated at 255 W (m^2)
n_s = ceil(A_c / A_panel); % No. of solar panels needed
Ps_rated = 255; % Rated max power of each module (W)
Ps_rated_total = ns * Ps_rated; % Total solar PV installed capacity in microgrid (W)
L_s = 25; % Average lifetime of PV array (y)

IC_s = (Ps_rated_total/1000) * 1800; % Assuming installed PV costs of $1800/kW
OM_s = (2/100)*IC_s; % Fixed annual O&M costs ($/y)

% Wind
h_hub = 50; % Desired hub height (m)
h_ref = 5; % Reference height at which hourly wind speed v_ref is measured (m)
eta_w = 0.25; % Wind turbine generator efficiency (obtained from manufacturer data)
A = 1521; % Wind turbine rotor swept area (m^2)
IC_w = (Pw_rated/1000) * 1500; % Capital costs depend on power rating - from IRENA
OM_w = (2/100)*IC_w; % Fixed O&M costs ($/y)
L_w = 20; % Average operating lifetime of turbine (y)

% Use historical hourly wind speed data over a year (or an average hourly profile for just one day)
v_ref = [];

%% CONVENTIONAL (FF) GENERATORS
% Microturbine (gas)

% Diesel generator 

%% STORAGE
% Upper and lower limits on battery state of charge depend on depth of discharge (DoD) and battery state of health (SoH)
% Need track of total no. of cycles to determine replacement at end of lifetime
% Final source used for all parameters: https://www.nrel.gov/docs/fy16osti/64987.pdf

% Lithium ion battery pack similar to Tesla Powerwall
% More recent source for Li-ion: https://www.tesla.com/powerwall
n_LI; % No. of Li-ion battery packs
cycles_LI; % Cumulative no. of cycles
E_init_LI = 13.5; % Initial rated max capacity (kWh)
nc_LI = 0; % Zero charge-discharge cycles at program start
SOC_min_LI = 10;
SOC_max_LI = 90; % DoD for LI = 80% (https://www.spiritenergy.co.uk/kb-batteries-understanding-batteries)
P_max_LI = 7000; % Max/peak rated power (W), 5 kW continuous
V_LI = 400; % Nominal cell voltage (in V)
IC1_LI = 520; % Initial capital capacity costs ($/kWh) 
L_LI = 15; % Lifetime assumed to be 5 y beyond warranty period (y)

% $209/kWh - Source: https://www.nrel.gov/docs/fy19osti/71714.pdf
%IC2_LI = 120; % Initial capital power costs ($/kW)
% OM_LI = 5; % O&M costs ($/kW-yr)
IC_LI = IC1_LI * E_init_LI; % Total initial capital costs for LI (upfront)
OM_LI = (2/100) * IC_LI; % Fixed O&M costs ($/y)

n_LA;
cycles_LA;
E_init_LA = 1.68; % Max rated capacity of each battery (kWh) 
SOC_min_LA = 40;
SOC_max_LA = 90; % DoD for LA = 50%
V_LA = 12; % Nominal cell voltage (in V) 
P_max_LA = 420; % Max rated power (W)
IC1_LA = 255; % ($/kWh)
IC2_LA = 125; % ($/kW) https://journals.uic.edu/ojs/index.php/JUR/article/view/7556
% OM_LA = 10; % Fixed O&M costs ($/kW-yr) https://journals.uic.edu/ojs/index.php/JUR/article/view/7556
IC_LA = IC1_LA * E_init_LA; % Total initial capital costs for LA
OM_LA = (2/100) * IC_LA; % Fixed O&M costs ($/y)
% Assume no variable O&M costs for WT/PV/BS - no. of cycles only affects % when replacement is needed

%% Constraints

% Overall system - 
% Power balance (1) Supply = demand
% Reliability (2) Upper limits on LLSP, LLF etc.
% (3) All numbers n_i >= 0

% Solar

% Wind

% Battery 
% (4) SOC_min <= SOC <= SOC_max
% (5) P_min <= P_b <= P_max/P_rated
% (6) Limits on battery temp too?

% MT and DE
% (7) P_DE and P_MT >= 0 (sanity check: need to enforce since these are set based on differences)
%% Objective functions 
% To be computed/recalculated inside the CVX program at each iteration/timestep

% Overall system lifecycle cost (annualized) 
% Net present cost
% GHG Emissions
% Reliability indices

%% Design optimization


%% Optimal dispatch of storage and generation
% Minimize overall operational costs and emissions of the system

