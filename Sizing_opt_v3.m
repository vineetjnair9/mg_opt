cvx_begin
    variable n_s integer
    variable n_w integer
    variable Eb_init 
    
    % Initialize variables
    expression P_s(8760); % (kW)
    expression P_w(8760); % (kW)
    expression P_DE(8760); % (kW)
    expression P_b_LI(8760); % (kW)
    expression P_RES(8760); % (kW)
    expression P_dump(8760); % Excess power (RES output) sent to dump loads/ground (kW)
    expression soc_LI(8761); % between 0 and 1
    expression DE_ON(8760); % Boolean variable to keep track of whether DE is online/committed
    % Lost load or Amount of demand response (curtailment/shifting) needed to balance MG
    expression P_lost(8760); 
    
    % LOAD ALL FIXED EXTERNAL/EXOGENOUS INPUT VARIABLES
    load('irradiance.mat'); % Solar insolation/irradiance (kW/m^2)
    load('temp.mat'); % Ambient air temperature (C)
    load('wind_speed.mat'); % (m/s)
    load('Load.mat'); % Hourly power (kW) demand profiles 
    % Hourly energy demand profile = power demand * time interval
    
    % COSTS
    i = 0.1; % Nominal interest/discount rate (from NREL LCOE explorer)
    f = 0.035; % Inflation/escalation rate (https://tradingeconomics.com/tanzania/inflation-cpi)
    r = (i-f)/(1+f); % Real discount/interest rate
    t_overall = 25; % Overall system lifetime (y) = lifetime of solar (longest surviving components)
    CRF = (r*(1+r)^t_overall)/((1+r)^t_overall - 1); % Capital recovery factor

    % Startup and shutdown costs of conventional generators
    DE_startup = 0;
    DE_shutdown = 0;

    % Variable O&M costs of DE
    Var_OM_DE = 0;

    % Fuel costs
    C_fuel_DE = 0; 
    CO2 = 0; 

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
    t_s = 25; % Lifetime of solar PV array (y)
    Ps_rated = 255; % Rated max power of each module (W)
    eta_pc = 1; % Power conditioning efficiency (= 1 or 100% since it's assumed that a perfect MPPT is used)
    Ps_rated_total = n_s * Ps_rated; % Total solar PV installed capacity in microgrid (W)
    IC_s = (Ps_rated_total/1000) * 1750; % Assuming installed PV costs of $1750/kW (source: IRENA costs report 2018, total installed costs)
    OM_s = (1/100)*IC_s; % Fixed annual O&M costs ($/y) - from Kaabeche et al. 2011a

    % WIND
    beta_w = 1/7; % Typical value of power law exponent for open land
    % Use data from this spec sheet: https://farm-energy.extension.org/wp-content/uploads/2019/04/3.5kW_Spec_Sheet.pdf
    h_ref = 2; % Installation height of TAHMO measurement stations (m)
    h_hub = 14.5; % Tower height to hub/nacelle (m)
    Pw_rated = 3.5; % Rated power of each WT (kW)
    Pw_rated_total = n_w * Pw_rated;
    IC_w = Pw_rated_total * 1500; % Assuming installed PV costs of $1500/kW (source: IRENA costs report 2018, total installed costs)
    OM_w = (3/100)*IC_w; % Fixed annual O&M costs ($/y)
    t_w = 20; % Lifetime of wind turbine system (y)
    v_c = 2.8; % Cut-in speed (m/s)
    v_r = 11; % Rated wind speed (m/s)
    v_f = 22; % Cut-out or failure speed (m/s)
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
    IC1_LI = 520; % Initial capital energy/capacity costs ($/kWh)
    IC2_LI = 120; % Initial capital power costs ($/kW) - source?
    IC_LI = IC1_LI * Eb_init + IC2_LI * P_max_LI; % Total initial capital costs for LI 
    OM_LI = (1/100) * IC_LI; % from Kaabeche et al. 2011a ($/y)
    RC_LI = 0.9*IC_LI; % Replacement cost at end of service period/lifetime
    Eb_init = Eb_init*3.6e6; % 1 W = 1 J/s
    E_C = Eb_init; % Actual capacity is initially at rated maximum (J)
    t_LI = 10; % Average lifetime of Li-ion battery (y) - Tesla powerwall warranty 

    % Diesel generator specs from Moshi et al. 2016
    P_DE_rated = 16; % Rated maximum power of diesel generator (kW) = approx. peak load * 1.5
    P_DE_min = 4.8; % Lower limit on power (kW) or k_gen * P_rated (k_gen = 0.3) from Fathima et al. 
    Ramp = (P_DE_rated*1000)/10; % Max ramp rate of DE (W/min) - NERC disturbance control standard 
    % Ramp_up =  0.006038 - 0.000003840*(P_DE_rated/10e3); % from ramp_rates.pdf (as % of nameplate capacity per minute)
    % Ramp_down = 0.006783 - 0.000004314*(P_DE_rated/10e3);
    IC_DE = 12500; % Initial capital cost ($)
    RC_DE = 11000; % Replacement cost ($)
    Fixed_OM_DE = (2/100)*IC_DE; % Fixed O&M costs ($/y)
    Variable_OM_DE = 0.24; % Variable O&M costs ($/h of operation online)
    SUC_DE = 0.45; % Start-up cost ($/kW)
    SDC_DE = 0.23; % Shut-down cost ($)
    t_DE = 5; % Average lifetime of diesel engine (y) - source?
    % t_DE = 15000; % Lifetime (h), from Moshi et al. 2016

    % Bidirectional inverter between AC and DC buses (can act as either rectifier or inverter)
    eta_inv = 0.9; % Inverter efficiency (Ogunjuyigbe et al. 2016)
    eta_rec = eta_inv; % Both rectifier and inverter assumed to have same parameters (Moshi et al. 2016)
    t_inv = 10; % Lifetime of converter 
    P_inv_rated = 16; % Maximum rated power of inverter (kW) - chosen to be same as that of backup generator
    IC_inv = 0.715 * P_inv_rated; % Assuming unit price of $0.715/W
    % Assume converter has zero maintenance costs

    % Start with delta_t = 1 h and then maybe discretize over smaller time intervals
    delta_t = 3600; % Calculation period (s) = 1 hour

    hour = 0; % Keeps track of the hour in the day

    % Normalize emissions wrt a base case where the DE is continuously run for the whole year at rated power - to meet all the load
    fuel_DE = 8760 * 0.06 * P_DE_rated + 0.0246 * sum(sum(Load));
    CO2 = 3.5 * fuel_DE;
    CO = 4.063288937 * sum(sum(Load));
    NOx = 18.85658039 * sum(sum(Load));
    SO2 = 0.007381439 * sum(sum(Load));
    VOC = 1.502443664 * sum(sum(Load));
    PM = 1.338208931 * sum(sum(Load));
    PM10 = 1.338208931 * sum(sum(Load));
    PM25 = 1.338208931 * sum(sum(Load));
    Emissions_base = CO2 + (CO + NOx + SO2 + VOC + PM + PM10 + PM25)/1000; % (kg)

    CO2 = 0;

    for t = 1:1:8760 % Simulate over one year with a time-step of 1 h

        Excess_demand = 0;

        % SOLAR - How does Ps_rated affect these calculations?
        I = irradiance(t); % Hourly solar irradiance (kW/m^2)
        T_a = temp(t);
        eta_PV = eta_r * (1 - 0.9 * beta_s * (I/I_PV_NOCT) * (T_c_NOCT - T_a_NOCT) - beta_s * (T_a - T_r)); % Tazvinga et al. 2013
        % Alternative model - Kaabeche et al. 2011b
        T_c = T_a + ((T_c_NOCT - 20)/0.8)*I; % Cell temp (C)
        % eta_PV = eta_r * eta_pc * (1 - beta_s*(T_c - T_r));
        % A_c = 1.625 * 1.019; % Surface area of 1 panel rated at 255 W (m^2)
        A_c = 120 * 78e-3 * 156e-3; % Surface area of 1 panel rated at 255 W (m^2)
        P_s(t) = n_s * eta_PV * A_c * I; % Solar power output (kW) 

    %     if (P_s(t) > Ps_rated_total)
    %         P_s(t) = Ps_rated_total;
    %     end

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
            P_w(t) = n_w * Pw_rated; % Wind power output (kW)
        end

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
        P_load = Load(day_number,hour); % (kW)

        if (P_RES(t) == P_load/eta_inv)
            continue
        elseif (P_RES(t) > P_load/eta_inv) % Excess supply/generation
            if (soc_LI(t) >= SOC_max_LI) % Battery can't be charged further
                P_dump(t) = P_RES(t) - P_load/eta_inv; % Send excess power to dump loads/ground
            else % Charge battery system
                if (P_RES(t) - P_load/eta_inv <= P_max_LI) 
                    P_b_LI(t) = -(P_RES(t) - P_load/eta_inv);
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
            Excess_demand = P_load/eta_inv - P_RES(t);
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

                if (Excess_demand/eta_rec <= P_DE_rated)
                    if (Excess_demand/eta_rec >= P_DE_min)
                        P_DE(t) = Excess_demand/eta_rec;
                        % No dumped load
                    else
                        P_DE(t) = P_DE_min;
                        P_dump(t) = P_DE_min - Excess_demand/eta_rec; 
                        % If we allow excess generator output to charge battery
                        if (soc_LI(t) < SOC_max_LI)
                            if (P_dump(t) <= P_max_LI)
                                if (cycled == 0) % Battery hasn't been cycled/operated yet before in this time step
                                    soc_LI(t+1) = soc_LI(t)*(1-delta_hour) + P_dump(t)*delta_t*(eta_overall/E_C); 
                                else % Battery has already been cycled so don't include self-discharge again
                                    soc_LI(t+1) = soc_LI(t+1) + P_dump(t)*delta_t*(eta_overall/E_C);
                                end
                                P_dump(t) = 0;
                            else
                                if (cycled == 0) % Battery hasn't been cycled/operated yet before in this time step
                                    soc_LI(t+1) = soc_LI(t)*(1-delta_hour) + P_max_LI*delta_t*(eta_overall/E_C); 
                                else % Battery has already been cycled so don't include self-discharge again
                                    soc_LI(t+1) = soc_LI(t+1) + P_max_LI*delta_t*(eta_overall/E_C);
                                end                           
                                P_dump(t) = P_dump(t) - P_max_LI;
                            end
                        end
                    end
                else
                    P_DE(t) = P_DE_rated;
                    P_lost(t) = Excess_demand/eta_rec - P_DE_rated;
                end

                if (t ~= 1 && DE_ON(t-1) == 0) % i.e. DE was offline in the last/most recent time step
                    DE_startup = DE_startup + SUC_DE;
                    % Ramp_actual = (Excess_Demand - P_DE(t-1))/60 % actual ramp rate (W/min)
                end
            else
                DE_ON(t) = 0; % Turn DE off
                if (t ~= 1 && DE_ON(t-1) == 1) % i.e. DE was online in the last/most recent time step
                    DE_shutdown = DE_shutdown + SDC_DE;
                end
            end
        end

        fuel_DE = 0.06 * P_DE_rated + 0.0246 * P_DE(t); % fuel consumption (L/h) where P is in kW - Kaabeche and Ibtiouen 2014
        C_fuel_DE = C_fuel_DE + 3.20 * (fuel_DE / 3.78541); % Diesel fuel cost assuming a price of $3.20/US liquid gallon (value for Tanzania from NREL ReOpt)
        % Alternatively, can model assume fuel consumption cost to be quadratic function of DE power output (Parisio et al. 2014 and others)
        % But need to find right coefficients by data fitting
        CO2 = CO2 + 3.5 * fuel_DE; % Total CO2 emissions (kg) using emission factor of 3.5 kg/L of diesel - source?

        cycled = 0; % Reset
    end

    if (DE_ON(1) == 1)
        DE_startup = DE_startup + SUC_DE;
    end

    % Upon completion of for loop i.e. after iterating through all hours of the year, C_fuel_DE now gives the total fuel cost over the 1st year of MG operation
    % It is assumed that on average, the annual fuel cost remains constant for all subsequent years of the system's lifetime
    Var_OM_DE = Variable_OM_DE * sum(DE_ON);

    % Total initial capital/installed cost
    IC = IC_s + IC_w + IC_LI + IC_DE + IC_inv;

    % Total annual recurring costs each year
    C_rec = OM_s + OM_w + OM_LI + Fixed_OM_DE + Var_OM_DE + C_fuel_DE;

    % Present worth of recurring costs 
    PW_rec = C_rec * (((1+f)/(1+i))*(((1+f)/(1+i))^t_overall - 1))/(((1+f)/(1+i)) - 1);

    % Battery replaced every 10 years
    i_adj_LI = (((1+i)^t_LI)/(1+f)^(t_LI - 1)) - 1; % Adjusted nominal interest rate
    PW_LI_rep = RC_LI * (((1+f)/(1+i_adj_LI))*(((1+f)/(1+i_adj_LI))^t_overall - 1))/(((1+f)/(1+i_adj_LI)) - 1);

    % DE replaced every 5 years
    i_adj_DE = (((1+i)^t_DE)/(1+f)^(t_DE - 1)) - 1; % Adjusted nominal interest rate
    PW_DE_rep = RC_DE * (((1+f)/(1+i_adj_DE))*(((1+f)/(1+i_adj_DE))^t_overall - 1))/(((1+f)/(1+i_adj_DE)) - 1);

    % Present worth of non-recurring costs 
    PW_nonrec = PW_LI_rep + PW_DE_rep;

    TNPC = IC + PW_rec + PW_nonrec; % Total (lifecycle) net present cost of system ($)
    TAC = TNPC * CRF; % Total annualized cost ($) 
    E_gen = sum(P_DE + P_RES); % Total annual electrical energy = P_gen * delta_t = P_gen (in kWh) since delta_t = 1 h
    LCOE = TAC/E_gen; % Energy cost/Levelized cost of electricity ($/kWh)

    % Base case for LCOE same as that for emissions - entire MG is run solely on a single DE 
    % Total annual recurring costs each year
    Var_OM_DE = Variable_OM_DE * 8760;
    fuel_DE = 8760 * 0.06 * P_DE_rated + 0.0246 * sum(sum(Load));
    C_fuel_DE = 3.20 * (fuel_DE / 3.78541);
    C_rec = Fixed_OM_DE + Var_OM_DE + C_fuel_DE;

    % Present worth of recurring costs 
    PW_rec = C_rec * (((1+f)/(1+i))*(((1+f)/(1+i))^t_overall - 1))/(((1+f)/(1+i)) - 1);

    % Present worth of non-recurring costs 
    PW_nonrec = PW_DE_rep;

    IC =  IC_DE + IC_inv;
    TNPC = IC + PW_rec + PW_nonrec; % Total (lifecycle) net present cost of system ($)
    TAC = TNPC * CRF; % Total annualized cost ($) 
    LCOE_base = TAC/sum(sum(Load)); % Energy cost/Levelized cost of electricity ($/kWh)

    % Emissions factors - EPA database (all in g)
    % Alternative source: Wu et al. 2014
    CO = 4.063288937 * sum(P_DE);
    NOx = 18.85658039 * sum(P_DE);
    SO2 = 0.007381439 * sum(P_DE);
    VOC = 1.502443664 * sum(P_DE);
    PM = 1.338208931 * sum(P_DE);
    PM10 = 1.338208931 * sum(P_DE);
    PM25 = 1.338208931 * sum(P_DE);
    Emissions = CO2 + (CO + NOx + SO2 + VOC + PM + PM10 + PM25)/1000;
    Emissions = Emissions/Emissions_base; % Always lies between 0 and 1

    % Reliability measure - % of load that's not not met
    DPSP = sum(P_lost)/sum(sum(Load));

    % Excess/dump energy as a fraction of total generation
    Dump = sum(P_dump)/E_gen;

    % Renewable penetration
    REF = sum(P_RES)/E_gen;

    % time = 1:1:8760;
    % plot(time, P_RES);

    Costs = [LCOE/LCOE_base Emissions DPSP Dump 1-REF];
    w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
    cost = Costs * w';
    
    minimize cost
    subject to
        n_s >= 1e-4, n_w >= 1e-4, Eb_init >= 1e-4
cvx_end