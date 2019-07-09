function [Pw] = Wind(v_ref, h_hub, h_ref, eta_w, A, v_rated, Pw_rated)
% v_ref is the measured hourly wind speed at the reference height
% Obtained from historical data for chosen project location
% Model from Tazvinga et al. 2014: Energy dispatch strategy for a photovoltaic–wind–diesel–battery hybrid power system

rho_air = 1.225; % Density of air (kg/m^3)
beta = 1/7; % Typical value of power law exponent for open land

% Account for variation of Cp with wind speed from manufacturer data sheets
% Source for wind turbine characteristics - https://www.enercon.de/fileadmin/Redakteur/Medien-Portal/broschueren/pdf/en/ENERCON_Produkt_en_06_2015.pdf
% Use values for 900 kW ENERCON E-44 system

% Power coefficient of wind turbine (from ENERCON data sheet)
C_p = [0;0;0.16;0.34;0.43;0.48;0.49;0.50;0.50;0.50;0.48;0.44;0.39;0.33;0.28;0.24;0.20;0.17;0.14;0.12;0.11;0.09;0.08;0.07;0.06];

v_hub = v_ref * (h_hub/h_ref)^beta; % Wind speed at hub calculated from measured speed at reference anemometer height
P_wind = 0.5 * eta_w * rho_air * Cp * A * v_hub^3;
t_w = 20; % Lifetime of wind turbine system (y)
v_cutin = 3; % Cut-in speed (m/s)
v_cutout = 31; % Cut-out speed (m/s)

% From Borhanazad et al. 2014: Optimization of micro-grid system using MOPSO
if (v_hub < v_cutin || v_hub > v_cutout)
    Pw = 0;
elseif (v_cutin <= v_hub && v_hub <= v_rated)
    Pw = P_wind;
    % Alternative formula
    % Pw = ((v_hub^3)*P_rated - P_rated*v_cutin^3)/(v_rated^3 - v_cutin^3);
else
    Pw = Pw_rated;
end


end
