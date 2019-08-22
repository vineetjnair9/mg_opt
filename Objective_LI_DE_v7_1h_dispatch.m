function [P_RES, P_s, P_w, P_b_LI, P_DE, Loads, P_dump, soc_LI] = Objective_LI_DE_v7_1h_dispatch(x)
% To compute actual dispatch used for optimal sizing

% Function performs optimal dispatch to output overall weighted cost function that needs to be minimized
% Variable inputs (that we minimize w.r.t) are: 
% no. of solar PV arrays, no. of wind turbines, initial battery capacity (in kWh)
n_s = x(1);
n_w = x(2);
Eb_init = x(3);

% Given fixed inputs: Wind speed, ambient external air temperature, solar radiation, electricity demand
% Have all of this data available over a whole year with hourly resolution

% LOAD ALL FIXED EXTERNAL/EXOGENOUS INPUT VARIABLES
load('irradiance.mat'); % Solar insolation/irradiance (kW/m^2)
load('temp.mat'); % Ambient air temperature (C)
load('wind_speed.mat'); % (m/s)
load('Load.mat'); % Hourly power (kW) demand profiles 
% Hourly energy demand profile = power demand * time interval

% Initialize variables
P_s = zeros(1,8760); % (kW) - DC
P_w = zeros(1,8760); % (kW) - AC
P_DE = zeros(1,8760); % (kW) - AC
P_b_LI = zeros(1,8760); % (kW) - DC
Loads = zeros(1,8760);
P_RES = zeros(1,8760); % (kW) - DC
P_dump = zeros(1,8760); % Excess power (RES output) sent to dump loads/ground as DC (kW)
soc_LI = zeros(1,8761); % between 0 and 1
DE_ON = zeros(1,8760); % Boolean variable to keep track of whether DE is online/committed during that hour
% Lost load or Amount of demand response (curtailment/shifting) needed to balance MG
P_lost = zeros(1,8760); % (kW) - AC

% COSTS
i = 0.09; % Nominal interest/discount rate (https://tradingeconomics.com/kenya/interest-rate)
f = 0.057; % Inflation/escalation rate (https://tradingeconomics.com/kenya/inflation-cpi)
r = (i-f)/(1+f); % Real discount/interest rate
t_overall = 25; % Overall system lifetime (y) = lifetime of solar (longest surviving components)
CRF = (r*(1+r)^t_overall)/((1+r)^t_overall - 1); % Capital recovery factor

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
Ps_rated = 255; % (W)
Ps_rated_total = n_s * Ps_rated;
beta_s = 0.0045; % Temp coefficient of module/generator efficiency (of Pmax)
t_s = 25; % Lifetime of solar PV array (y)
eta_pc = 1; % Power conditioning efficiency (= 1 or 100% since it's assumed that a perfect MPPT is used)
IC_s = (Ps_rated_total/1000) * 1210; % Assuming installed PV costs of $1210/kW (source: IRENA costs report 2018, total installed costs)
% likely to be a slight underestimate since this is the global weighted
% average for larger utility scale projects
OM_s = (1/100)*IC_s; % Fixed annual O&M costs ($/y) - from Kaabeche et al. 2011a

% WIND
beta_w = 1/7; % Typical value of power law exponent for open land
% Use data from this spec sheet: https://farm-energy.extension.org/wp-content/uploads/2019/04/3.5kW_Spec_Sheet.pdf
h_ref = 1; % Installation height of TAHMO measurement stations (m) - wind speed measured approx. at sea level
h_hub = 14.5; % Tower height to hub/nacelle (m)
Pw_rated = 3.5; % (kW)
Pw_rated_total = n_w * Pw_rated;
IC_w = Pw_rated_total * 1500; % Global weighted average total installed costs of onshore wind ($/kW) - IRENA 2018

OM_w = (3/100)*IC_w; % Fixed annual O&M costs ($/y)
t_w = 20; % Lifetime of wind turbine system (y)
v_c = 2.8; % Cut-in or start-up wind speed (m/s)
v_r = 11; % Rated wind speed (m/s)
v_f = 22; % Cut-out or failure/braking speed (m/s)
A = (1/(v_c^2 - v_r^2)) * (v_c*(v_c + v_r) - 4*v_c*v_r*((v_c + v_r)/(2*v_r))^3);
B = (1/(v_c^2 - v_r^2)) * (4*(v_c + v_r)*((v_c + v_r)/(2*v_r))^3 - 3*(v_c + v_r));
C = (1/(v_c^2 - v_r^2)) * (2 - 4*((v_c + v_r)/(2*v_r))^3);

% Li-ion battery system - from https://www.nrel.gov/docs/fy16osti/64987.pdf 
% Updated parameters for powerwall - https://www.tesla.com/en_GB/powerwall
SOC_min_LI = 0.1;
SOC_max_LI = 0.9; % DoD for LI = 80% (https://www.spiritenergy.co.uk/kb-batteries-understanding-batteries) 
delta = 0.075; % Self-discharge rate of battery (7.5%/month)
P_max_LI = 3.68; % Max continuous real power (kW), 5 kW peak
delta_hour = delta/30.5; % Hourly self-discharge rate (%/h)
eta_overall = 0.90; % Battery round trip coulombic efficiency (accounts for both charge & discharge) - from Hu et al. 2017
cycles_LI = 0; % Keeps track of cumulative no. of charge/discharge cycles so far
cycled = 0; % Keeps track of whether the BS has already been self-discharged in the current time step
% IC1_LI = 300; % Initial energy/capacity capital costs ($/kWh) - NREL
% IC2_LI = 120; % Initial power capital costs ($/kW) 
% IC_LI = IC1_LI * Eb_init + IC2_LI * P_max_LI; % Total initial capital costs for LI 
IC_LI = 300 * Eb_init; % Assuming a price of $300/kWh
OM_LI = (1/100) * IC_LI; % from Kaabeche et al. 2011a ($/y)
RC_LI = 0.9*IC_LI; % Replacement cost at end of service period/lifetime - assumed to be 90% of IC since new installation and permitting fees aren't incurred
Eb_init = Eb_init*3.6e6; % 1 W = 1 J/s
E_C = Eb_init; % Actual capacity is initially at rated maximum (J)
%t_LI = 15; % Average lifetime of Li-ion battery (y) - 5 years beyond Tesla powerwall warranty period
t_LI = 5475; % cycles

% Diesel generator specs from Moshi et al. 2016
P_DE_rated = 16; % Rated maximum power of diesel generator (kW) = approx. peak load * 1.5
P_DE_min = 4.8; % Lower limit on power (kW) or k_gen * P_rated (k_gen = 0.3) from Fathima et al. 
Ramp = (P_DE_rated*1000)/10; % Max ramp rate of DE (W/min) - NERC disturbance control standard, DG must be able to reach rated capacity within 10 min
% Ramp_up =  0.006038 - 0.000003840*(P_DE_rated/10e3); % from ramp_rates.pdf (as % of nameplate capacity per minute)
% Ramp_down = 0.006783 - 0.000004314*(P_DE_rated/10e3);
IC_DE = 12500; % Initial capital cost for 16 kW rated DE ($) - slightly higher than replacement 
RC_DE = 11000; % Replacement cost ($)
Fixed_OM_DE = (2/100)*IC_DE; % Fixed O&M costs ($/y)
Variable_OM_DE = 0.24; % Variable O&M costs ($/h of operation online)
SUC_DE = 0.45; % Start-up cost ($)
SDC_DE = 0.23; % Shut-down cost ($)
% t_DE = 14; % Average lifetime of diesel engine (y) - avg expected life expectancy (https://www.wpowerproducts.com/news/diesel-engine-life-expectancy/)
t_DE = 15000; % Lifetime (h), from Moshi et al. 2016

% Bidirectional converter between AC and DC buses (can act as either rectifier or inverter)
eta_inv = 0.9; % Inverter efficiency (Ogunjuyigbe et al. 2016)
eta_rec = eta_inv; % Both rectifier and inverter assumed to have same parameters (Moshi et al. 2016)
t_inv = 20; % Lifetime of converter (y) - Moshi et al. 2016
P_inv_rated = 12.52; % Maximum rated power of inverter (kW)
IC_inv = 2800; % ($)
% Assume converter has zero maintenance costs

% Start with delta_t = 1 h and then maybe discretize over smaller time intervals
delta_t = 3600; % Calculation period (s) = 1 hour

hour = 0; % Keeps track of the hour in the day

for t = 1:1:8760 % Simulate over one year with a time-step of 1 h
    
    Excess_demand = 0;
    
    % SOLAR
    I = irradiance(t); % Hourly solar irradiance (kW/m^2)
    T_a = temp(t);
    eta_PV = eta_r * (1 - 0.9 * beta_s * (I/I_PV_NOCT) * (T_c_NOCT - T_a_NOCT) - beta_s * (T_a - T_r)); % Tazvinga et al. 2013
    % Alternative model - Kaabeche et al. 2011b
    % T_c = T_a + ((T_c_NOCT - 20)/0.8)*I; % Cell temp (C)
    % eta_PV = eta_r * eta_pc * (1 - beta_s*(T_c - T_r));
    A_c = 120 * 78e-3 * 156e-3; % Surface area of 1 panel rated at 255 W (m^2)
    P_s(t) = n_s * eta_PV * A_c * I; % Solar power output (kW) - DC
    
    if (P_s(t) > Ps_rated_total)
        P_s(t) = Ps_rated_total;
    end
    
    % WIND
    v_ref = wind_speed(t);
    v_hub = v_ref * (h_hub/h_ref)^beta_w; % Wind speed at hub calculated from measured speed at reference anemometer height
    
    if (v_hub < v_c || v_hub > v_f)
        P_w(t) = 0;
    elseif (v_c <= v_hub && v_hub <= v_r)
        P_w(t) = Pw_rated_total * (A + B * v_hub + C * v_hub^2); % Source?
        % Alternative formula - Borhanazad et al. 2014
        % Pw = ((v_hub^3)*P_rated - P_rated*v_cutin^3)/(v_rated^3 - v_cutin^3);
        % Also include commonly known formula using Cp and swept area? 
        % But need to know variation of Cp with TSR for my specific WT model      
    else
        P_w(t) = Pw_rated_total; % Wind power output (kW) - AC
    end
    
    % Rectify AC power from WT to DC
    P_RES(t) = P_w(t)*eta_rec + P_s(t); % Total renewable power output available at DC bus (kW)
    
    % BATTERY
    % Assume battery starts out fully charged (at max safe SOC)
    if (t == 1)
        soc_LI(t) = 0.9;
    end
    
    if (soc_LI(t) > SOC_max_LI)
        soc_LI(t) = SOC_max_LI;
    end
    
    if (soc_LI(t) < SOC_min_LI)
        soc_LI(t) = SOC_min_LI;
    end
    
    if (t ~= 1 && soc_LI(t) == SOC_max_LI && soc_LI(t-1) < SOC_max_LI) % Complete 1 charge-discharge cycle
        cycles_LI = cycles_LI + 1; % Increment cycle number
        E_C = Eb_init*(1 - (cycles_LI*0.055e-3)); % Capacity fading with cycling - https://www.nrel.gov/docs/fy16osti/64987.pdf
    end
      
    day_number = ceil(t/24);
    if (hour >= 24)
        hour = 0;
    end
    hour = hour + 1;
    P_load = Load(day_number,hour); % (kW) - AC
    Loads(t) = Load(day_number,hour);
       
    if (P_RES(t) == P_load/eta_inv) % Since RES power needs to inverted to meet AC load
        continue
    elseif (P_RES(t) > P_load/eta_inv) % Excess supply/generation
        if (soc_LI(t) >= SOC_max_LI) % Battery can't be charged further
            P_dump(t) = P_RES(t) - P_load/eta_inv; % Send excess power to dump loads/ground - DC
        else % Charge battery system
            if (P_RES(t) - P_load/eta_inv <= P_max_LI) 
                P_b_LI(t) = -(P_RES(t) - P_load/eta_inv); % (kW)
                soc_LI(t+1) = soc_LI(t)*(1-delta_hour) - P_b_LI(t)*1000*delta_t*(eta_overall/E_C);
                cycled = 1;
            else
                P_b_LI(t) = -P_max_LI;
                soc_LI(t+1) = soc_LI(t)*(1-delta_hour) - P_b_LI(t)*1000*delta_t*(eta_overall/E_C);
                P_dump(t) = P_RES(t) - P_load/eta_inv - P_max_LI;
                cycled = 1;
            end
        end
    else % Insufficient supply/generation
        Excess_demand = P_load/eta_inv - P_RES(t); % DC
        if (soc_LI(t) > SOC_min_LI) % Discharge battery
            if (Excess_demand <= P_max_LI)
                P_b_LI(t) = Excess_demand;
                soc_LI(t+1) = soc_LI(t)*(1-delta_hour) - P_b_LI(t)*1000*delta_t*(eta_overall/E_C);
                Excess_demand = 0; % All load has been met/satisfied
                cycled = 1;
            else
                P_b_LI(t) = P_max_LI;
                soc_LI(t+1) = soc_LI(t)*(1-delta_hour) - P_b_LI(t)*1000*delta_t*(eta_overall/E_C);
                Excess_demand = Excess_demand - P_max_LI;
                cycled = 1;
            end
        elseif (soc_LI(t) <= SOC_min_LI || Excess_demand > 0) % Need to use backup DE to meet excess load
            DE_ON(t) = 1; % Turn DE on
            % Power produced by DE is already AC - so can directly meet load
            Excess_demand = Excess_demand * eta_inv; % AC
            if (Excess_demand <= P_DE_rated)
                if (Excess_demand >= P_DE_min)
                    P_DE(t) = Excess_demand;
                    % No dumped load
                else
                    P_DE(t) = P_DE_min;
                    P_dump(t) = (P_DE_min - Excess_demand) * eta_rec; % (kW) - DC
                    % If we allow excess generator output to charge battery
                    % DE output must 1st be rectified to DC before charging battery
                    if (soc_LI(t) < SOC_max_LI)
                        if (P_dump(t) <= P_max_LI)
                            if (cycled == 0) % Battery hasn't been cycled/operated yet before in this time step
                                soc_LI(t+1) = soc_LI(t)*(1-delta_hour) + P_dump(t)*1000*delta_t*(eta_overall/E_C); 
                            else % Battery has already been cycled so don't include self-discharge again
                                soc_LI(t+1) = soc_LI(t+1) + P_dump(t)*1000*delta_t*(eta_overall/E_C);
                            end
                            P_dump(t) = 0;
                        else
                            if (cycled == 0) % Battery hasn't been cycled/operated yet before in this time step
                                soc_LI(t+1) = soc_LI(t)*(1-delta_hour) + P_max_LI*1000*delta_t*(eta_overall/E_C); 
                            else % Battery has already been cycled so don't include self-discharge again
                                soc_LI(t+1) = soc_LI(t+1) + P_max_LI*1000*delta_t*(eta_overall/E_C);
                            end                           
                            P_dump(t) = P_dump(t) - P_max_LI;
                        end
                    end
                end
            else
                P_DE(t) = P_DE_rated;
                P_lost(t) = Excess_demand - P_DE_rated;
            end

        else
            DE_ON(t) = 0; % Turn DE off
        end
    end
    
    cycled = 0; % Reset
end

end