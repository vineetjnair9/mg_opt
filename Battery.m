function [soc, E_new] = Battery(soc_prev, Pb, E_init, delta_t)
% Generalized battery model
% TODO: Find more detailed/realistic models specific to the battery type used in final design (Pb-acid, Li-ion, NiCd, NiMH etc.)

% P_batt is the power 'output' from the battery
% E_init = Initial rated capacity/size of the battery storage (BS) system = Design variable (in Ah) 
% TODO: Using dummy values for now, need to find correct battery parameters

%% Simpler model obtained from Tazvinga et al. 2014: Energy dispatch strategy for a photovoltaic–wind–diesel–battery hybrid power system
% Assumes the max battery capacity to be constant and always = rated capacity 
% So calculating absolute rather than relative SoC
eta_overall = 1; % Battery round trip efficiency (accounts for both charge & discharge)
alpha = eta_overall/E_init;
soc = soc_prev - Pb*alpha;

%% Model obtained from Wu et al. 2014: Dynamic economic dispatch of a microgrid: Mathematical models and solution algorithm
% 1st attempt: Keep max battery capacity constant
% TODO: Later refine model to account for capacity fading (reducing) with increased cycling and use
% TODO: Will need to keep track of cumulative number of cycles and update E_C at each iteration (charge/discharge cycle)

eta_c = 1; % charging efficiency
eta_d = 1; % discharging efficiency
Pb_max = 1; % Max power output limit of the battery

% E_init = Total initial capacity of battery system (in Ah)
E_C = E_init ; % Total current capacity of battery during calculation period delta_t

delta = 0.01; % Self-discharge rate of battery (%/h)

if (Pb >= 0) % Charging regime
    soc = (1-delta)*soc_prev + (Pb*delta_t*eta_c)/E_C;
    E_new = E_C + (Pb/Pb_max)*E_C; % Rem
else 
    soc = (1-delta)*soc_prev - (Pb*delta_t)/(E_C*eta_d);
    E_new = E_C - (Pb/Pb_max)*E_C;
end

% TODO: Need to incorporate degradation costs of ES components 

end
