% Optimal dispatch of DE and BS on a day-ahead basis

load('irradiance.mat'); % Solar insolation/irradiance (kW/m^2)
load('temp.mat'); % Ambient air temperature (C)
load('wind_speed.mat'); % (m/s)
load('Load.mat'); % Hourly power (kW) demand profiles

% Need some form of receding horizon MPC since in reality, controller
% doesn't have access to future climate/load data!

cvx_begin
    variable P_DE
    variable P_b 
    
    for t = 1:1:23
        
        
        
    
    
    