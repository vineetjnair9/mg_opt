%% Overall optimization problem
% Define the global variables
delta_t = 1; % time step or calculation period (s)

% Solar
A_c = 1; % Total surface area of PV array (m^2) 
eta_PV = 0.2; % Overall efficiency of PV array system
% Latitude and longitude of microgrid location as user inputs

% Wind
h_hub = 30; % Desired hub height (m)
h_ref = 5; % Reference height at which hourly wind speed v_ref is measured (m)
eta_w = 0.25; % Wind turbine generator efficiency (obtained from manufacturer data)
C_p = 0.3; % Power coefficient of wind turbine 
A = 100; % Wind turbine rotor swept area (m^2)

%% Cost functions

%% Constraints
% Upper and lower limits on battery state of charge 
SOC_min = 20;
SOC_max = 100;

%% Design optimization
% Size capacities of generation and storage components to minimize net present costs (NPC) 
% as well as annualized costs over the system's lifecycle