function [cost, soc_LI, DPSP, P_dump, P_lost, P_load, P_w, P_s, P_RES] = Dispatch_obj_LI_DE_v4_DPSP(P_DE, P_b_LI)
% Optimal dispatch of DE and BS on a day-ahead basis - Using inbuilt MATLAB solvers
% Initial attempt: No MPC, static optimization over finite horizon based on certain/known historical load/climate data i.e. similar to LQR
% But in reality, need some form of receding horizon MPC since controller doesn't have access to future climate/load data!

% Assuming that any excess DE power cannot go towards charging battery
% i.e. purpose of DE is only to serve as backup to meet excess load
% P_DE and P_b_LI are 1x24 vectors

load('irradiance.mat'); % Solar insolation/irradiance (kW/m^2)
load('temp.mat'); % Ambient air temperature (C)
load('wind_speed.mat'); % (m/s)
load('load_peaky.mat'); % Hourly power (kW) demand profile for a single representative day

% Optimal sizing results using 5 equal weights - from v2
% Oversize RES (round up n_s and n_w to nearest integer) and storage capacity        
% Using DG = 16 kW
% n_s = 15.1187;
% n_w = 4.3061;
% Eb_init = 106.5333; % (kWh)

% Using DG = 8 kW
n_s = 8.0095;
n_w = 4.6226;
Eb_init = 118.1596; % (kWh)

% Discretize 1 day into 1 h time intervals 
delta_t = 3600; % Calculation period (s) = 1 h

P_s = zeros(1,24); % (kW)
P_w = zeros(1,24); % (kW)
P_RES = zeros(1,24); % (kW)
P_dump = zeros(1,24); % Excess power (RES output) sent to dump loads/ground (kW)
soc_LI = zeros(1,25); 
DE_ON = zeros(1,24); % Boolean variable to keep track of whether DE is online/committed
P_lost = zeros(1,24); % Lost load AC
P_load = zeros(1,24);

% Startup and shutdown costs of conventional generators
DE_startup = 0;
DE_shutdown = 0;

% Fuel costs
C_fuel_DE = 0; 

% Fixed system parameters
% SOLAR 
% PV panel specs from https://www.mitsubishielectricsolar.com/images/uploads/documents/specs/MLU_spec_sheet_250W_255W.pdf
% For 255 W system - Monocrystalline Si
eta_r = 0.154; % Manufacturer rated efficiency
T_r = 25; % Standard test conditions (reference cell temperature)
I_PV_NOCT = 0.8; % Hourly solar irradiance at NOCT - Normal/nominal operating cell temp (in kW/m^2) 
T_c_NOCT = 45.7; % PV cell temp at NOCT (Celsius)
T_a_NOCT = 20; % Ambient temp at NOCT (Celsius)
beta_s = 0.0045; % Temp coefficient of module/generator efficiency (midway between the Pmax and Voc coefficients)
Ps_rated = 255; % Rated max power of each module (W)
eta_pc = 1; % Power conditioning efficiency (= 1 or 100% since it's assumed that a perfect MPPT is used)
Ps_rated_total = n_s * Ps_rated; % Total solar PV installed capacity in microgrid (W)
IC_s = (Ps_rated_total/1000) * 1210; % Assuming installed PV costs of $1210/kW (source: IRENA costs report 2018, total installed costs)
OM_s = (1/100)*(IC_s/365); % Fixed daily O&M costs ($/y) - from Kaabeche et al. 2011a

% WIND
beta_w = 1/7; % Typical value of power law exponent for open land
% Use data from this spec sheet: https://farm-energy.extension.org/wp-content/uploads/2019/04/3.5kW_Spec_Sheet.pdf
h_ref = 1; % Installation height of TAHMO measurement stations (m)
h_hub = 14.5; % Tower height to hub/nacelle (m)
Pw_rated = 3.5; % Rated power of each WT (kW)
Pw_rated_total = n_w * Pw_rated;
IC_w = Pw_rated_total * 1500; % Global weighted average total installed costs of onshore wind ($/kW) - IRENA 2018
OM_w = (3/100)*(IC_w/365); % Fixed daily O&M costs ($/y)
v_c = 2.8; % Cut-in speed (m/s)
v_r = 11; % Rated wind speed (m/s)
v_f = 22; % Cut-out or failure speed (m/s)
A = (1/(v_c^2 - v_r^2)) * (v_c*(v_c + v_r) - 4*v_c*v_r*((v_c + v_r)/(2*v_r))^3);
B = (1/(v_c^2 - v_r^2)) * (4*(v_c + v_r)*((v_c + v_r)/(2*v_r))^3 - 3*(v_c + v_r));
C = (1/(v_c^2 - v_r^2)) * (2 - 4*((v_c + v_r)/(2*v_r))^3);

% Li-ion battery system - from https://www.nrel.gov/docs/fy16osti/64987.pdf 
SOC_min_LI = 0.1;
SOC_max_LI = 0.9; % DoD for LI = 80% (https://www.spiritenergy.co.uk/kb-batteries-understanding-batteries) 
delta = 0.075; % Self-discharge rate of battery (%/month)
P_max_LI = 3.68; % Max continuous real power (kW), 5 kW peak
IC_LI = 300 * Eb_init;
OM_LI = (1/100) * (IC_LI/365); % Fixed daily O&M cost from Kaabeche et al. 2011a ($/y)
delta_hour = delta/30.5; % Hourly self-discharge rate (%/h)
eta_overall = 0.90; % Battery round trip coulombic efficiency (accounts for both charge & discharge) - from Hu et al. 2017
cycles_LI = 0; % Keeps track of cumulative no. of charge/discharge cycles so far
Eb_init = Eb_init*3.6e6; % 1 W = 1 J/s
E_C = Eb_init; % Actual capacity is initially at rated maximum (J)

% Diesel generator specs from Moshi et al. 2016
P_DE_rated = 8; % Rated maximum power of diesel generator (kW) = approx. peak load * 1.5
P_DE_min = 0.3*P_DE_rated; % Lower limit on power (kW) or k_gen * P_rated (k_gen = 0.3) from Fathima et al. 
Ramp = (P_DE_rated*1000)/10; % Max ramp rate of DE (W/min) - NERC disturbance control standard 
% Ramp_up =  0.006038 - 0.000003840*(P_DE_rated/10e3); % from ramp_rates.pdf (as % of nameplate capacity per minute)
% Ramp_down = 0.006783 - 0.000004314*(P_DE_rated/10e3);
IC_DE = 12500; % Initial capital cost ($)
Fixed_OM_DE = (2/100)*(IC_DE/365); % Fixed daily O&M costs ($/d)
Variable_OM_DE = 0.24; % Variable O&M costs ($/h of operation online)
SUC_DE = 0.45; % Start-up cost ($/kW)
SDC_DE = 0.23; % Shut-down cost ($)

% Bidirectional inverter between AC and DC buses (can act as either rectifier or inverter)
eta_inv = 0.9; % Inverter efficiency (Ogunjuyigbe et al. 2016)
eta_rec = eta_inv; % Both rectifier and inverter assumed to have same parameters (Moshi et al. 2016)
t_inv = 20; % Lifetime of converter 
P_inv_rated = 16; % Maximum rated power of inverter (kW) - chosen to be same as that of backup generator

P_load = load_peaky;

% Normalize emissions wrt a base case where the DE is continuously run for the whole day - to meet all the load
CO2 = 0.649 * sum(P_load);
CO = 4.063288937 * sum(P_load);
NOx = 18.85658039 * sum(P_load);
SO2 = 0.007381439 * sum(P_load);
VOC = 1.502443664 * sum(P_load);
PM = 1.338208931 * sum(P_load);
PM10 = 1.338208931 * sum(P_load);
PM25 = 1.338208931 * sum(P_load);
Emissions_base = CO2 + (CO + NOx + SO2 + VOC + PM + PM10 + PM25)/1000; % (kg)
  
CO2 = 0;
    
for t = 1:1:24
    % SOLAR
    I = irradiance(t); % Hourly solar irradiance (kW/m^2)
    T_a = temp(t);
    eta_PV = eta_r * (1 - 0.9 * beta_s * (I/I_PV_NOCT) * (T_c_NOCT - T_a_NOCT) - beta_s * (T_a - T_r)); % Tazvinga et al. 2013
    A_c = 120 * 78e-3 * 156e-3; % Surface area of 1 panel rated at 255 W (m^2)
    P_s(t) = n_s * eta_PV * A_c * I; % Solar power output (kW) 

    if (P_s(t) > Ps_rated_total)
        P_s(t) = Ps_rated_total;
    end
    
    % WIND
    v_ref = wind_speed(t);
    v_hub = v_ref * (h_hub/h_ref)^beta_w; % Wind speed at hub calculated from measured speed at reference anemometer height

    % From Borhanazad et al. 2014: Optimization of micro-grid system using MOPSO
    if (v_hub < v_c || v_hub > v_f)
        P_w(t) = 0;
    elseif (v_c <= v_hub && v_hub <= v_r)
        P_w(t) = n_w * Pw_rated * (A + B * v_hub + C * v_hub^2);
    else
        P_w(t) = n_w * Pw_rated; % Wind power output (kW)
    end

    P_RES(t) = P_w(t)*eta_rec + P_s(t); % Total renewable power output available at DC bus (kW)

    if (t == 1)
        soc_LI(t) = 0.9;
    end

    if (t ~= 1 && soc_LI(t) == SOC_max_LI && soc_LI(t-1) < SOC_max_LI) % Complete 1 charge-discharge cycle
        cycles_LI = cycles_LI + 1; % Increment cycle number
        E_C = Eb_init*(1 - (cycles_LI*0.055e-3)); % Capacity fading with cycling - https://www.nrel.gov/docs/fy16osti/64987.pdf
    end

    % Lower limit on P_DE or not
    if (P_DE(t) >= P_DE_min) % DE online during hour t (P_DE_min or 1e-4)
        DE_ON(t) = 1;
    end
    
    % To enforce lower limit constraint on P_DE
    if (P_DE(t) < P_DE_min)
        P_DE(t) = 0;
    end

    if (DE_ON(1) == 1)
        DE_startup = DE_startup + SUC_DE;
    end

    if (t ~= 1 && DE_ON(t) == 1 && DE_ON(t-1) == 0) % i.e. DE was offline in the last/most recent time step
        DE_startup = DE_startup + SUC_DE;  
    end

    if (t ~= 1 && DE_ON(t-1) == 1) % i.e. DE was online in the last/most recent time step
        DE_shutdown = DE_shutdown + SDC_DE;
    end

    fuel_DE = 0.08145 * P_DE_rated + 0.246 * P_DE(t);
    C_fuel_DE = C_fuel_DE + 3.20 * (fuel_DE / 3.78541);
    CO2 = CO2 + 0.649 * P_DE(t);

    % Battery dynamics
    % P_b_LI +ve -> Discharging
    % P_b_LI -ve -> Charging
    soc_LI(t+1) = soc_LI(t)*(1-delta_hour) - P_b_LI(t)*1000*delta_t*(eta_overall/E_C);

    if (P_RES(t) == P_load(t)/eta_inv)
        continue
    elseif (P_RES(t) > P_load(t)/eta_inv) % Excess supply
        % Need to rectify AC DE output before sending to DC bus for dump loads
        P_dump(t) = max(0, P_RES(t) + P_b_LI(t) + P_DE(t)*eta_rec - P_load(t)/eta_inv); % Sign of P_b_LI takes care of both charge/discharge case
    else % Excess demand (AC)
        Excess_demand = P_load(t) - (P_RES(t) * eta_inv); 
        P_lost(t) = max(0, Excess_demand - P_DE(t) - (P_b_LI(t) * eta_inv)); % AC
        P_dump(t) = max(0, (P_DE_min - Excess_demand + (P_b_LI(t) * eta_inv)) * eta_rec); % DC
    end
end

CO = 4.063288937 * sum(P_DE);
NOx = 18.85658039 * sum(P_DE);
SO2 = 0.007381439 * sum(P_DE);
VOC = 1.502443664 * sum(P_DE);
PM = 1.338208931 * sum(P_DE);
PM10 = 1.338208931 * sum(P_DE);
PM25 = 1.338208931 * sum(P_DE);
Emissions = CO2 + (CO + NOx + SO2 + VOC + PM + PM10 + PM25)/1000;
Emissions = Emissions/Emissions_base; % Always lies between 0 and 1

Var_OM_DE = Variable_OM_DE * sum(DE_ON);
Cost = C_fuel_DE + Var_OM_DE + OM_s + OM_w + OM_LI + Fixed_OM_DE + DE_startup + DE_shutdown;
E_gen = sum(P_DE * eta_rec + P_RES);
COE = Cost/sum(P_load); % Unit cost of electricity ($/kWh) - considering only real-time, operational costs

% Base case for LCOE same as that for emissions - entire MG is run solely on a single DE 
Var_OM_DE = Variable_OM_DE * 24;
fuel_DE = 24 * 0.08145 * P_DE_rated + 0.246 * sum(P_load);
C_fuel_DE = 3.20 * (fuel_DE / 3.78541);
DE_startup = SUC_DE;
DE_shutdown = SDC_DE;
Cost_base = C_fuel_DE + Var_OM_DE + Fixed_OM_DE + DE_startup + DE_shutdown;
COE_base = Cost_base/sum(P_load);

DPSP = sum(P_lost)/sum(P_load);
Dump = sum(P_dump)/E_gen;
REF = sum(P_RES)/E_gen;

% No longer need to min DPSP since it's a constraint
Costs = [COE/COE_base Emissions DPSP Dump 1-REF];
w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = Costs * w';

end
  