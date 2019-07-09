function [Ps] = Solar(A_c, I_PV)
% Simple model to calculate ohourly energy output from PV array
% From Tazvinga et al. 2014: Energy dispatch strategy for a photovoltaic–wind–diesel–battery hybrid power system
% A_c = Total collector area of the PV array (m^2) 
% eta_PV = Electrical conversion efficiency of the array
% I_PV = Hourly solar radiation (kWh/m^2) - incident on the PV array per hour

% Solar cell efficiency and power output also dependent on  external environmental conditions 
% Use ambient temperature data along with actual horizontal irradiation / surface irradiance data for chosen location in TZ

% PV panel specs from https://www.mitsubishielectricsolar.com/images/uploads/documents/specs/MLU_spec_sheet_250W_255W.pdf
% For 255 W system - Monocrystalline Si
% Formula for varying solar PV efficiency from Hove and Tazvinga et al. 2012 
% Also used Tazvinga et al. 2013

eta_r = 0.154; % Manufacturer rated efficiency
T_r = 25; % Standard test conditions (room temperature in Celsius)
I_PV_NOCT = 0.8; % Hourly solar radiation at NOCT - Normal operating cell temp (in kWh/m^2)
T_C_NOCT = 45.7; % PV cell temp at NOCT (Celsius)
T_a_NOCT = 20; % Ambient temp at NOCT
beta = -0.45e-2; % Temp coefficient of module efficiency (midway between the Pmax and Voc coefficients)
eta_PV = eta_r * (1 - 0.9 * (I_PV/I_PV_NOCT) * (T_C_NOCT - T_a_NOCT) - beta * (T_a - T_r));

Ps = eta_PV * A_c * I_PV;

end