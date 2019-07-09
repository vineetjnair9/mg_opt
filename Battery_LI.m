function soc = Battery_LI(soc_prev, Pb, delta_t, cycles_LI)
% Generalized battery model
% P_b is the total power 'output' from the battery storage system
% E_init = Initial rated max capacity/size of the battery storage (BS) system = Design variable (in Ah) 
Pb_max = 7000; % Max rated power output limit of the battery (W)

%% Simpler model obtained from Tazvinga et al. 2014: Energy dispatch strategy for a photovoltaic–wind–diesel–battery hybrid power system
% Assumes the max battery capacity to be constant and always = rated capacity 
% So calculating absolute rather than relative SoC

eta_overall = 0.91; % Battery round trip coulombic efficiency (accounts for both charge & discharge)
% from Hu et al. 2017
E_init = 13500; % Wh
delta = 0.075; % Self-discharge rate of battery (%/month)
E_C = E_init*(1- (cycles_LI*0.0055)/100); % Capacity fading with cycling - https://www.nrel.gov/docs/fy16osti/64987.pdf
soc = soc_prev*(1-delta) - Pb*delta_t*(eta_overall/E_C);
% E_new = E_C - (Pb/Pb_max)*E_C;

%% Model obtained from Wu et al. 2014: Dynamic economic dispatch of a microgrid: Mathematical models and solution algorithm
% 1st attempt: Keep max battery capacity constant
% TODO: Later refine model to account for capacity fading (reducing) with increased cycling and use
% TODO: Will need to keep track of cumulative number of cycles and update E_C at each iteration (charge/discharge cycle)

% Most parameters from https://batteryuniversity.com/learn/article/secondary_batteries

% Assume eta_c = eta_d = sqrt(eta_overall)
eta_c = 0.954; % charging efficiency
eta_d = 0.954; % discharging efficiency

% E_init = Total initial capacity of battery system (in Ah)
E_C = E_init ; % Total current capacity of battery during calculation period delta_t

delta = 0.075; % Self-discharge rate of battery (%/month) -- need conversion depending on delta_t units

if (Pb < 0) % Charging regime
    soc = (1-delta)*soc_prev - (Pb*delta_t*eta_c)/E_C;
    E_new = E_C - (Pb/Pb_max)*E_C; % Remaining energy
else % Discharge regime (+ve power output from battery)
    soc = (1-delta)*soc_prev - (Pb*delta_t)/(E_C*eta_d);
    E_new = E_C - (Pb/Pb_max)*E_C;
end

end

