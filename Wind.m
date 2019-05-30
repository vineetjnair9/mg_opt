function [Pw] = Wind(v_ref, h_hub, h_ref, eta_w, C_p, A, v_rated, P_rated, v_cutin, v_cutout, x_w)
% v_ref is the measured hourly wind speed at the reference height
% Obtained from historical data for chosen project location
% Model from Tazvinga et al. 2014: Energy dispatch strategy for a photovoltaic–wind–diesel–battery hybrid power system
% x_w = Number of wind turbines turned on (committed) in the time interval being considered

rho_air = 1.225; % Density of air (kg/m^3)
beta = 1/7; % Typical value of power law exponent for open land
v_hub = v_ref * (h_hub/h_ref)^beta;
P_wind = 0.5 * eta_w * rho_air * C_p * A * v_hub^3;

% TODO: Account for power-limiting below cut-in and cut-out speed
% From Borhanazad et al. 2014: Optimization of micro-grid system using MOPSO
if (v_hub < v_cutin || v_hub > v_cutout)
    Pw = P_wind;
elseif (v_cutin <= v_hub && v_hub <= v_rated)
    Pw = P_wind;
    % Alternative formula
    % Pw = ((v_hub^3)*P_rated - P_rated*v_cutin^3)/(v_rated^3 - v_cutin^3);
else
    Pw = P_rated;
end
end
