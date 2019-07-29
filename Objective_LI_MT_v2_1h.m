function [cost] = Objective_LI_MT_v2_1h(x)
% With modified strategy for excess demand scenario/case (elseif statement & not self-discharging twice in the same time-step)
% Normalized LCOE and emissions to always lie between 0 and 1

n_s = x(1);
n_w = x(2);
Eb_init = x(3);

% LOAD ALL FIXED EXTERNAL/EXOGENOUS INPUT VARIABLES
load('irradiance.mat'); % Solar insolation/irradiance (kW/m^2)
load('temp.mat'); % Ambient air temperature (C)
load('wind_speed.mat'); % (m/s)
load('Load.mat'); % Hourly power (kW) demand profiles 

% Initialize variables
P_s = zeros(1,8760); % (kW) - DC
P_w = zeros(1,8760); % (kW) - AC
P_MT = zeros(1,8760); % (kW) - AC
P_b_LI = zeros(1,8760); % (kW) - DC
P_RES = zeros(1,8760); % (kW) - DC
P_dump = zeros(1,8760); % Excess power (RES output) sent to dump loads/ground as DC (kW)
soc_LI = zeros(1,8761); % between 0 and 1
MT_ON = zeros(1,8760); % Boolean variable to keep track of whether MT is online/committed during that hour
P_lost = zeros(1,8760); % (kW) - AC

% COSTS
i = 0.09; % Nominal interest/discount rate (https://tradingeconomics.com/kenya/interest-rate)
f = 0.057; % Inflation/escalation rate (https://tradingeconomics.com/kenya/inflation-cpi)
r = (i-f)/(1+f); % Real discount/interest rate
t_overall = 25; % Overall system lifetime (y) = lifetime of solar (longest surviving components)
CRF = (r*(1+r)^t_overall)/((1+r)^t_overall - 1); % Capital recovery factor

% Startup and shutdown costs of conventional generators
MT_startup = 0;
MT_shutdown = 0;

% Fuel costs
C_fuel_MT = 0; 

% Fixed system parameters
% SOLAR 
eta_r = 0.154; % Manufacturer rated efficiency
T_r = 25; % Standard test conditions (reference cell temperature)
I_PV_NOCT = 0.8; % Hourly solar irradiance at NOCT - Normal/nominal operating cell temp (in kW/m^2) 
T_c_NOCT = 45.7; % PV cell temp at NOCT (Celsius)
T_a_NOCT = 20; % Ambient temp at NOCT (Celsius)
beta_s = 0.0045; % Temp coefficient of module/generator efficiency (of Pmax)
t_s = 25; % Lifetime of solar PV array (y)
Ps_rated = 255; % Rated max power of each module (W)
eta_pc = 1; % Power conditioning efficiency (= 1 or 100% since it's assumed that a perfect MPPT is used)
Ps_rated_total = n_s * Ps_rated; % Total solar PV installed capacity in microgrid (W)
IC_s = (Ps_rated_total/1000) * 1210; % Assuming installed PV costs of $1210/kW (source: IRENA costs report 2018, total installed costs)
OM_s = (1/100)*IC_s; % Fixed annual O&M costs ($/y) - from Kaabeche et al. 2011a

% WIND
beta_w = 1/7; % Typical value of power law exponent for open land
% Use data from this spec sheet: https://farm-energy.extension.org/wp-content/uploads/2019/04/3.5kW_Spec_Sheet.pdf
h_ref = 1; % Installation height of TAHMO measurement stations (m) - wind speed measured approx. at sea level
h_hub = 14.5; % Tower height to hub/nacelle (m)
Pw_rated = 3.5; % Rated power of each WT (kW)
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
IC_LI = 300 * Eb_init; % Assuming a price of $300/kWh
OM_LI = (1/100) * IC_LI; % from Kaabeche et al. 2011a ($/y)
RC_LI = 0.9*IC_LI; % Replacement cost at end of service period/lifetime - assumed to be 90% of IC since new installation and permitting fees aren't incurred
Eb_init = Eb_init*3.6e6; % 1 W = 1 J/s
E_C = Eb_init; % Actual capacity is initially at rated maximum (J)
t_LI = 15; % Average lifetime of Li-ion battery (y) - 5 years beyond Tesla powerwall warranty period

% MT specs from https://www.energy.gov/sites/prod/files/2016/09/f33/CHP-Microturbines_0.pdf
P_MT_rated = 16; % Rated maximum power of MT (kW) = approx. peak load * 1.5
P_MT_min = 4.8; % Lower limit on power (kW) or k_gen * P_rated (k_gen = 0.3) from Fathima et al. 
Ramp = (P_MT_rated*1000)/10; % Max ramp rate of MT (W/min) - NERC disturbance control standard, DG must be able to reach rated capacity within 10 min
IC_MT = 3220 * P_MT_rated; % Initial capital cost - https://www.energy.gov/sites/prod/files/2016/09/f33/CHP-Microturbines_0.pdf
RC_MT = 0.9*IC_MT; % Replacement cost ($)
Fixed_OM_MT = (2/100)*IC_MT; % Fixed O&M costs ($/y)
Variable_OM_MT = 0.013; % Variable O&M costs ($/kWh)
C_NG = 2.19; % Cost of natural gas in $/MMBtu - https://markets.businessinsider.com/commodities/natural-gas-price
SUC_MT = 0.45; % Start-up cost ($) - assum similar startup/shutdown costs to MT
SDC_MT = 0.23; % Shut-down cost ($)
t_MT = 40000; % Lifetime (h), from https://www.epa.gov/sites/production/files/2015-07/documents/catalog_of_chp_technologies_section_5._characterization_-_microturbines.pdf

% Bidirectional converter between AC and DC buses (can act as either rectifier or inverter)
eta_inv = 0.9; % Inverter efficiency (Ogunjuyigbe et al. 2016)
eta_rec = eta_inv; % Both rectifier and inverter assumed to have same parameters (Moshi et al. 2016)
t_inv = 20; % Lifetime of converter (y) - Moshi et al. 2016
P_inv_rated = 16; % Maximum rated power of inverter (kW) - chosen to be larger than peak load (1.5x)
IC_inv = 2 * P_inv_rated * 1000; % Assuming unit price of $2/W - https://www.nrel.gov/docs/fy19osti/72399.pdf

% Start with delta_t = 1 h and then maybe discretize over smaller time intervals
delta_t = 3600; % Calculation period (s) = 1 hour

hour = 0; % Keeps track of the hour in the day

% Normalize emissions wrt a base case where the MT is continuously run for the whole year at rated power - to meet all the load
CO2 = 0.631 * sum(sum(Load));
CO = 2.851087583 * sum(P_MT);
NOx = 20.88408858 * sum(P_MT);
SO2 = 0.003009766 * sum(P_MT);
VOC = 0.604000601  * sum(P_MT);
PM = 0.046974355 * sum(P_MT);
PM10 = 0.046974355 * sum(P_MT);
PM25 = 0.046974355 * sum(P_MT);
Emissions_base = CO2 + (CO + NOx + SO2 + VOC + PM + PM10 + PM25)/1000; % (kg)

CO2 = 0;

for t = 1:1:8760 % Simulate over one year with a time-step of 1 h
    
    Excess_demand = 0;
    
    % SOLAR
    I = irradiance(t); % Hourly solar irradiance (kW/m^2)
    T_a = temp(t);
    eta_PV = eta_r * (1 - 0.9 * beta_s * (I/I_PV_NOCT) * (T_c_NOCT - T_a_NOCT) - beta_s * (T_a - T_r)); % Tazvinga et al. 2013
    A_c = 120 * 78e-3 * 156e-3; % Surface area of 1 panel rated at 255 W (m^2)
    P_s(t) = n_s * eta_PV * A_c * I; % Solar power output (kW) - DC
    
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
        % Alternative formula
        % Pw = ((v_hub^3)*P_rated - P_rated*v_cutin^3)/(v_rated^3 - v_cutin^3);
    else
        P_w(t) = n_w * Pw_rated; % Wind power output (kW) - AC
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
        elseif (soc_LI(t) <= SOC_min_LI || Excess_demand > 0) % Need to use backup MT to meet excess load
            MT_ON(t) = 1; % Turn MT on
            % Power produced by MT is already AC - so can directly meet load
            Excess_demand = Excess_demand * eta_inv; % AC
            if (Excess_demand <= P_MT_rated)
                if (Excess_demand >= P_MT_min)
                    P_MT(t) = Excess_demand;
                    % No dumped load
                else
                    P_MT(t) = P_MT_min;
                    P_dump(t) = (P_MT_min - Excess_demand) * eta_rec; % (kW) - DC
                    % If we allow excess generator output to charge battery
                    % MT output must 1st be rectified to DC before charging battery
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
                P_MT(t) = P_MT_rated;
                P_lost(t) = Excess_demand - P_MT_rated;
            end
                        
            if (t ~= 1 && MT_ON(t-1) == 0) % i.e. MT was offline in the last/most recent time step
                MT_startup = MT_startup + SUC_MT;
                % Ramp_actual = (Excess_Demand - P_MT(t-1))/60 % actual ramp rate (W/min)
            end
        else
            MT_ON(t) = 0; % Turn MT off
            if (t ~= 1 && MT_ON(t-1) == 1) % i.e. MT was online in the last/most recent time step
                MT_shutdown = MT_shutdown + SDC_MT;
            end
        end
    end
    
    % Also need to account for (varying) MT efficiency while calculating its fuel consumption!    
    % fuel_MT = 0.06 * P_MT_rated + 0.0246 * P_MT(t); % fuel consumption (L/h) where P is in kW - Kaabeche and Ibtiouen 2014
    % C_fuel_MT = C_fuel_MT + 3.20 * (fuel_MT / 3.78541); % Diesel fuel cost assuming a price of $3.20/US liquid gallon (value for east Africa from NREL ReOpt)
    
    CO2 = CO2 + 0.631 * P_MT(t); % Usign emission coefficient of 0.631 kg/kWh - https://www.energy.gov/sites/prod/files/2016/09/f33/CHP-Microturbines_0.pdf
    cycled = 0; % Reset
end

if (MT_ON(1) == 1)
    MT_startup = MT_startup + SUC_MT;
end

Var_OM_MT = Variable_OM_MT * sum(P_MT);
fuel_MT = (0.84/61) * P_MT_rated * sum(MT_ON); % NG fuel consumption (MMBtu) 
% assuming 0.84 MMTBtu/h of operation for a 61 kW MT
C_fuel_MT = fuel_MT * C_NG;

% Total initial capital/installed cost
IC = IC_s + IC_w + IC_LI + IC_MT + IC_inv;

% Total annual recurring costs each year
C_rec = OM_s + OM_w + OM_LI + Fixed_OM_MT + Var_OM_MT + C_fuel_MT + MT_startup + MT_shutdown;

% Present worth of recurring costs 
PW_rec = C_rec * (((1+f)/(1+i))*(((1+f)/(1+i))^t_overall - 1))/(((1+f)/(1+i)) - 1);

% Battery replaced every 1400 cycles
t_rep_LI = cycles_LI/t_LI; % No. of years between battery replacements
i_adj_LI = (((1+i)^t_rep_LI)/(1+f)^(t_rep_LI - 1)) - 1; % Adjusted nominal interest rate
PW_LI_rep = RC_LI * (((1+f)/(1+i_adj_LI))*(((1+f)/(1+i_adj_LI))^t_overall - 1))/(((1+f)/(1+i_adj_LI)) - 1);
    
% Determine frequency of MT replacement based on actual no. of operating hours per year
t_MT_rep = (t_MT/sum(MT_ON)); % No. of years after which MT needs to be replaced
i_adj_MT = (((1+i)^t_MT_rep)/(1+f)^(t_MT_rep - 1)) - 1; % Adjusted nominal interest rate
PW_MT_rep = RC_MT * (((1+f)/(1+i_adj_MT))*(((1+f)/(1+i_adj_MT))^t_overall - 1))/(((1+f)/(1+i_adj_MT)) - 1);

% Present worth of non-recurring costs 
PW_nonrec = PW_LI_rep + PW_MT_rep;

TNPC = IC + PW_rec + PW_nonrec; % Total (lifecycle) net present cost of system ($)
TAC = TNPC * CRF; % Total annualized cost ($) 
LCOE = TAC/sum(sum(Load)); % Energy cost/Levelized cost of electricity ($/kWh) = total annual costs / total load served

% Base case for LCOE same as that for emissions - entire MG is run solely on a single MT
% Total annual recurring costs each yeaR
Var_OM_MT = Variable_OM_MT * 8760;
C_fuel_MT = fuel_MT * (0.84/61) * P_MT_rated * 8760;
MT_startup = SUC_MT;
MT_shutdown = SDC_MT;
C_rec = Fixed_OM_mt + Var_OM_MT + C_fuel_MT + MT_startup + MT_shutdown;

% Present worth of recurring costs 
PW_rec = C_rec * (((1+f)/(1+i))*(((1+f)/(1+i))^t_overall - 1))/(((1+f)/(1+i)) - 1);

t_MT_rep = (t_MT/8760); % No. of years after which MT needs to be replaced
i_adj_MT = (((1+i)^t_MT_rep)/(1+f)^(t_MT_rep - 1)) - 1; % Adjusted nominal interest rate
PW_MT_rep = RC_MT * (((1+f)/(1+i_adj_MT))*(((1+f)/(1+i_adj_MT))^t_overall - 1))/(((1+f)/(1+i_adj_MT)) - 1);

% Present worth of non-recurring costs 
PW_nonrec = PW_MT_rep;

IC =  IC_MT + IC_inv;
TNPC = IC + PW_rec + PW_nonrec; % Total (lifecycle) net present cost of system ($)
TAC = TNPC * CRF; % Total annualized cost ($) 
LCOE_base = TAC/sum(sum(Load)); % Energy cost/Levelized cost of electricity ($/kWh)

% Emissions factors - EPA database (all in g)
CO = 2.851087583 * sum(P_MT);
NOx = 20.88408858 * sum(P_MT);
SO2 = 0.003009766 * sum(P_MT);
VOC = 0.604000601  * sum(P_MT);
PM = 0.046974355 * sum(P_MT);
PM10 = 0.046974355 * sum(P_MT);
PM25 = 0.046974355 * sum(P_MT);
Emissions = CO2 + (CO + NOx + SO2 + VOC + PM + PM10 + PM25)/1000; % Total (kg)

% Reliability measure - % of load that's not not met
DPSP = sum(P_lost)/sum(sum(Load));

% Excess/dump energy as a fraction of total generation
E_gen = sum(P_MT * eta_rec + P_RES); % Total annual DC electrical energy output = P_gen * delta_t = P_gen (in kWh) since delta_t = 1 h
Dump = sum(P_dump)/E_gen;

% Renewable penetration
REF = sum(P_RES)/E_gen;

Costs = [LCOE/LCOE_base Emissions/Emissions_base DPSP Dump 1-REF];
w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = Costs * w';

end