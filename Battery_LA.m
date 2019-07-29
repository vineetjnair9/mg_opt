function [soc, E_new] = Battery_LA(soc_prev, P_b, E_init, delta_t, cycles_LA)
% Generalized battery model
% TODO: Find more detailed/realistic models specific to the battery type used in final design (Pb-acid, Li-ion, NiCd, NiMH etc.)

% P_batt is the total power 'output' from the battery storage system
% E_init = Initial rated capacity/size of the battery storage (BS) system = Design variable (in Ah) 
Pb_max = 420; % Max rated power output limit of the battery (W) - Charging/discharging

%% Simpler model obtained from Tazvinga et al. 2014: Energy dispatch strategy for a photovoltaic–wind–diesel–battery hybrid power system

eta_overall = 0.8; % Battery round trip efficiency (accounts for both charge & discharge)
delta = 0.05; % Self-discharge rate of battery (%/month)
E_C = E_init*(1- (cycles_LA*0.0214)/100); % Capacity fading with cycling - https://www.nrel.gov/docs/fy16osti/64987.pdf
soc = soc_prev*(1-delta) - P_b*delta_t*(eta_overall/E_C);
% E_new = E_C - (Pb/Pb_max)*E_C;

%% Model obtained from Wu et al. 2014: Dynamic economic dispatch of a microgrid: Mathematical models and solution algorithm
% 1st attempt: Keep max battery capacity constant
% TODO: Later refine model to account for capacity fading (reducing) with increased cycling and use
% TODO: update E_C at each iteration (charge/discharge cycle) according to total/cumulative no. of cycles thus far  

% eta_c = 0.8944; % charging efficiency
% eta_d = 0.8944; % discharging efficiency
% 
% % E_init = Initial maximum storage capacity (in Ah)
% E_C = E_init ; % Total current capacity of battery during calculation period delta_t
% 
% delta = 0.05; % Self-discharge rate of battery (%/h)
% 
% if (Pb < 0) % Charging regime
%     soc = (1-delta)*soc_prev - (Pb*delta_t*eta_c)/E_C;
%     E_new = E_C - (Pb/Pb_max)*E_C; % Remaining energy
% else % Discharge regime (+ve power output from battery)
%     soc = (1-delta)*soc_prev - (Pb*delta_t)/(E_C*eta_d);
%     E_new = E_C - (Pb/Pb_max)*E_C;
% end

end
