%% Daily dispatch power plot - w/ v2

P_DE = sol.P_DE;
P_b_LI = sol.P_b_LI;

[cost, soc_LI, DPSP, P_dump, P_lost, P_load, P_w, P_s, P_RES] = Dispatch_obj_LI_DE_v3(P_DE', P_b_LI');

P_DE_rated = 8; % (kW)
P_max_LI = 3.68; % (kW)

lw = 1.25;

t = 1:1:24;
figure(1)
plot(t,P_RES,'LineWidth',lw);
hold on;
plot(t, P_s,'LineWidth',lw);
plot(t, P_w,'LineWidth',lw);
plot(t, P_b_LI,'LineWidth',lw);
plot(t, P_DE,'LineWidth',lw);
plot(t, P_load,'LineWidth',lw);
plot(t, P_lost,'--','LineWidth',lw);
plot(t, P_dump,'--','LineWidth',lw);
% plot(t, P_max_LI.*ones(1,24),'c:','LineWidth',0.6);
% plot(t, -P_max_LI.*ones(1,24),'c:','LineWidth',0.6);
% plot(t, P_DE_rated.*ones(1,24),'m:','LineWidth',0.6);
set(gca,'FontSize',18);

xlabel('Time of the day (in h)');
ylabel('Power (in kW)');
legend({'P_{RES}','P_s','P_w','P_{LI}','P_{DE}','P_{load}','P_{lost}','P_{dump}'},'Orientation','horizontal');
xlim([1 24]);
%ylim([-2 25]);

figure(2)
plot(t, soc_LI(1:24).*100,'LineWidth',lw);
xlabel('Time of the day (in h)');
ylabel('SOC of LI BS (in %)')
set(gca,'FontSize',18);
xlim([1 24]);

%% Vary SOC plot

[Costs, cost, soc_LI, DPSP, P_dump, P_lost, P_load, P_w, P_s, P_RES] = Dispatch_obj_LI_DE_v5_varySOC(P_DE_soc(9,:), P_LI_soc(9,:), 0.9);

P_DE_rated = 8; % (kW)
P_max_LI = 3.68; % (kW)

lw = 1.25;

t = 1:1:72;
figure(1)
hold on;
plot(t,P_RES,'LineWidth',lw);
hold on;
plot(t, P_s,'LineWidth',lw);
plot(t, P_w,'LineWidth',lw);
plot(t, P_LI_soc(9,:),'LineWidth',lw);
plot(t, P_DE_soc(9,:),'LineWidth',lw);
plot(t, P_load,'LineWidth',lw);
plot(t, P_lost,'--','LineWidth',lw);
plot(t, P_dump,'--','LineWidth',lw);
% plot(t, P_max_LI.*ones(1,24),'c:','LineWidth',0.6);
% plot(t, -P_max_LI.*ones(1,24),'c:','LineWidth',0.6);
% plot(t, P_DE_rated.*ones(1,24),'m:','LineWidth',0.6);
set(gca,'FontSize',18);

xlabel('Time of the day (in h)');
ylabel('Power (in kW)');
legend({'P_{RES}','P_s','P_w','P_{LI}','P_{DE}','P_{load}','P_{lost}','P_{dump}'},'Orientation','horizontal');
xlim([1 72]);

figure(2)
plot(t, soc_LI(1:72).*100,'LineWidth',lw);
xlabel('Time of the day (in h)');
ylabel('SOC of LI BS (in %)')
set(gca,'FontSize',18);
xlim([1 72]);

%% Actual dispatch with opt size - for whole year
% Too noisy - may need to do some kind of smoothing/averaging

x = [15.1187    4.3061  106.5333];
[P_RES, P_s, P_w, P_b_LI, P_DE, Loads, P_dump, soc_LI] = Objective_LI_DE_v7_1h_dispatch(x);

time = 1:1:8760;

figure(1)
plot(time, P_RES);
hold on
plot(time, P_s);
plot(time, P_w);
plot(time, P_b_LI);
plot(time, P_DE);
plot(time, Loads);
%plot(t, P_lost);
plot(time, P_dump);
hold off
xlabel('Time of the year (in h)');
ylabel('Power (in kW)')
legend('P_{RES}','P_s','P_w','P_{LI}','P_{DE}','P_{load}','P_{dump}')
set(gca,'FontSize',18);

figure(2)
plot(time, soc_LI(1:8760).*100);
xlabel('Time of the year (in h)');
ylabel('SOC of LI BS (in %)')
set(gca,'FontSize',18);

%% Actual dispatch with opt size - for one day / week / month

x = [15.1187    4.3061  106.5333];
[P_RES, P_s, P_w, P_b_LI, P_DE, Loads, P_dump, soc_LI] = Objective_LI_DE_v7_1h_dispatch(x);

P_DE_rated = 16; % (kW)
P_max_LI = 3.68; % (kW)

t_end = 24;
time = 1:1:t_end;

lw = 1.25;

figure(1)
plot(time, P_RES(1:t_end),'LineWidth',lw);
hold on
plot(time, P_s(1:t_end),'LineWidth',lw);
plot(time, P_w(1:t_end),'LineWidth',lw);
plot(time, P_b_LI(1:t_end),'LineWidth',lw);
plot(time, P_DE(1:t_end),'LineWidth',lw);
plot(time, Loads(1:t_end),'--','LineWidth',lw);
plot(time, P_dump(1:t_end),'--','LineWidth',lw);
plot(time, P_max_LI.*ones(1,24),'c:','LineWidth',lw);
plot(time, -P_max_LI.*ones(1,24),'c:','LineWidth',lw);
plot(time, P_DE_rated.*ones(1,24),'m:','LineWidth',lw);
hold off
xlabel('Time of the day (in h)');
ylabel('Power (in kW)');
legend({'P_{RES}','P_s','P_w','P_{LI}','P_{DE}','P_{load}','P_{dump}'},'Orientation','horizontal');
set(gca,'FontSize',18);
%set(gca,'DefaultLineLineWidth',5)
xlim([1 t_end]);

figure(2)
plot(time, soc_LI(1:t_end).*100,'LineWidth',lw);
xlabel('Time of the month (in h)');
ylabel('SOC of LI BS (in %)');
xlim([1 t_end]);
ylim([0 100]);
set(gca,'FontSize',18);
%set(gca,'DefaultLineLineWidth',2)

%% Sizing varying parameter plots
% Varying DE rating
params = 0:1:20; % in kW
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = Costs * w';

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',DPSP,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
legend('Normalized emissions','DPSP','Dump','Weighted cost')
xlabel('Backup DE rated power (kW)');
ylabel('Objective value [-]');

figure(2)
hold on
plot(params',REF,'.-');
plot(params',LCOE,'.-');
legend('REF','Normalized LCOE');
xlabel('Backup DE rated power (kW)');
ylabel('Objective value [-]');

figure(3)
plot(params',x(:,3),'.-');
xlabel('Backup DE rated power (kW)');
ylabel('BS capacity (kWh)');

figure(4)
hold on
plot(params',DE_hours,'.-');
plot(params',BS_cycles,'.-');
xlabel('Backup DE rated power (kW)');
legend('No. of online DE hours per year [h]','No. of BS cycles per year');

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(5)
hold on
plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Backup DE rated power (kW)');
ylabel('Total rated power (kW)');
legend('PV capacity','WT capacity','Total RES capacity')

%% Varying DE fuel price
params = linspace(0.2,12,20); % in $/gal
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = Costs * w';

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
legend('Normalized emissions','Dump','Weighted cost')
xlabel('Diesel fuel price [$/gal]');
ylabel('Objective value [-]');

figure(2)
hold on
plot(params',REF,'.-');
plot(params',LCOE,'.-');
legend('REF','Normalized LCOE');
xlabel('Diesel fuel price [$/gal]');
ylabel('Objective value [-]');

figure(3)
plot(params',x(:,3),'.-');
xlabel('Diesel fuel price [$/gal]');
ylabel('BS capacity (kWh)');

figure(4)
hold on
plot(params',DE_hours,'.-');
plot(params',BS_cycles,'.-');
xlabel('Diesel fuel price [$/gal]');
legend('No. of online DE hours per year [h]','No. of BS cycles per year');

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(5)
hold on
plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Diesel fuel price [$/gal]');
ylabel('Total rated power (kW)');
legend('PV capacity','WT capacity','Total RES capacity')

%% Varying interest rate
params = linspace(0,0.2,20); % in (%)
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = Costs * w';

params = params.*100;

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
legend('Normalized emissions','Dump','Weighted cost')
xlabel('Nominal interest rate [%]');
ylabel('Objective value [-]');

figure(2)
hold on
plot(params',REF,'.-');
plot(params',LCOE,'.-');
legend('REF','Normalized LCOE');
xlabel('Nominal interest rate [%]');
ylabel('Objective value [-]');

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(3)
hold on
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Nominal interest rate [%]');
ylabel('Total rated power (kW)');
legend('WT capacity','Total RES capacity')

%% Varying inflation rate
params = linspace(0,0.15,20); 
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = Costs * w';

params = params.*100;

% LCOE values on a very different scale so need to plot separately
% figure(1)
% hold on
% plot(params',Emissions,'.-');
% plot(params',Dump,'.-');
% plot(params',cost,'.-');
% legend('Normalized emissions','Dump','Weighted cost')
% xlabel('Inflation rate (%)');
% ylabel('Objective value [-]');
% 
% figure(2)
% hold on
% plot(params',REF,'.-');
% plot(params',LCOE,'.-');
% legend('REF','Normalized LCOE');
% xlabel('Inflation rate (%)');
% ylabel('Objective value [-]');
% 
% figure(3)
% plot(params',x(:,3),'.-');
% xlabel('Inflation rate (%)');
% ylabel('BS capacity (kWh)');

figure(4)
hold on
plot(params',DE_hours,'.-');
%plot(params',BS_cycles,'.-');
xlabel('Inflation rate (%)');
ylabel('No. of online DE hours per year [h]');
%legend('No. of online DE hours per year [h]','No. of BS cycles per year');

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(5)
hold on
%plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Inflation rate (%)');
ylabel('Total rated power (kW)');
legend('WT capacity','Total RES capacity')

%% Varying BS price
params = linspace(50,300,20); 
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = Costs * w';

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
legend('Normalized emissions','Dump','Weighted cost')
xlabel('BS cost [$/kWh]');
ylabel('Objective value [-]');

figure(2)
hold on
plot(params',REF,'.-');
plot(params',LCOE,'.-');
legend('REF','Normalized LCOE');
xlabel('BS cost [$/kWh]');
ylabel('Objective value [-]');

figure(3)
plot(params',x(:,3),'.-');
xlabel('BS cost [$/kWh]');
ylabel('BS capacity (kWh)');

figure(4)
hold on
plot(params',DE_hours,'.-');
plot(params',BS_cycles,'.-');
xlabel('BS cost [$/kWh]');
legend('No. of online DE hours per year [h]','No. of BS cycles per year');

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(5)
hold on
plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('BS cost [$/kWh]');
ylabel('Total rated power (kW)');
legend('PV capacity','WT capacity','Total RES capacity')

%% Varying weight w1
params = linspace(0,1,11); 
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);
cost = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
    rem = (1-param)/4;
    w = [param rem rem rem rem];
    cost(i) = Costs(i,:) * w';
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
plot(params',LCOE,'.-');
legend('Normalized emissions','Dump','Weighted cost','Normalized LCOE')
xlabel('Value of weight w_1');
ylabel('Objective value [-]');

figure(2)
yyaxis left
plot(params',x(:,3),'.-');
xlabel('Value of weight w_1');
ylabel('BS capacity (kWh)');
yyaxis right
plot(params',BS_cycles,'.-');
ylabel('No. of BS cycles per year');
legend('BS capacity','Annual BS cycles');

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(3)
hold on
plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Value of weight w_1');
ylabel('Total rated power (kW)');
legend('PV capacity','WT capacity','Total RES capacity')

%% Varying weight w2
params = linspace(0,1,11); 
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);
cost = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
    rem = (1-param)/4;
    w = [rem param rem rem rem];
    cost(i) = Costs(i,:) * w';
end
LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
plot(params',REF,'.-');
plot(params',LCOE,'.-');
legend('Normalized emissions','Dump','Weighted cost','REF','Normalized LCOE')
xlabel('Value of weight w_2');
ylabel('Objective value [-]');

figure(2)
yyaxis left
plot(params',x(:,3),'.-');
xlabel('Value of weight w_2');
ylabel('BS capacity (kWh)');
yyaxis right
plot(params',DE_hours,'.-');
ylabel('No. of online DE hours per year [h]')
legend('BS capacity', 'DE hours')

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(3)
hold on
plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Value of weight w_2');
ylabel('Total rated power (kW)');
legend('PV capacity','WT capacity','Total RES capacity')

%% Varying w4
params = linspace(0,1,11); 
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);
cost = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
    rem = (1-param)/4;
    w = [rem rem rem param rem];
    cost(i) = Costs(i,:) * w';
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
plot(params',REF,'.-');
plot(params',LCOE,'.-');
legend('Normalized emissions','Dump','Weighted cost','REF','Normalized LCOE')
xlabel('Value of weight w_4');
ylabel('Objective value [-]');

figure(2)
yyaxis left
plot(params',x(:,3),'.-');
xlabel('Value of weight w_4');
ylabel('BS capacity (kWh)');
yyaxis right
plot(params',DE_hours,'.-');
ylabel('No. of online DE hours per year [h]')
legend('BS capacity', 'DE hours')

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(3)
hold on    
plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Value of weight w_4');
ylabel('Total rated power (kW)');
legend('PV capacity','WT capacity','Total RES capacity')

%% Varying w3
params = linspace(0,1,11); 
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);
cost = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
    rem = (1-param)/4;
    w = [rem rem param rem rem];
    cost(i) = Costs(i,:) * w';
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',cost,'.-');
plot(params',DPSP,'.-');
legend('Normalized emissions','Weighted cost','DPSP')
xlabel('Value of weight w_3');
ylabel('Objective value [-]');

figure(2)
hold on
plot(params',REF,'.-');
plot(params',LCOE,'.-');
plot(params',Dump,'.-');
legend('REF','Normalized LCOE','Dump');
xlabel('Value of weight w_3');
ylabel('Objective value [-]');

figure(3)
yyaxis left
plot(params',x(:,3),'.-');
xlabel('Value of weight w_3');
ylabel('BS capacity (kWh)');
yyaxis right
plot(params',DE_hours,'.-');
ylabel('No. of online DE hours per year [h]')
legend('BS capacity', 'DE hours')

RES_rated = x(:,1).*0.255 + x(:,2).*3.5;
figure(4)
hold on
plot(params',x(:,1).*0.255,'.-');
plot(params',x(:,2).*3.5,'.-');
plot(params',RES_rated,'.-');
xlabel('Value of weight w_3');
ylabel('Total rated power (kW)');
legend('PV capacity','WT capacity','Total RES capacity')

%% Varying w5
params = linspace(0,1,11); 
dim = size(params,2);
Costs = zeros(dim,5);
DE_hours = zeros(dim,1);
BS_cycles = zeros(dim,1);
cost = zeros(dim,1);

% Recording objective values for optimal solutions
for i = 1:1:dim
    param = params(i);
    [Costs(i,:),DE_hours(i),BS_cycles(i)] = Objective_LI_DE_v10_1h_ind_VaryParams_output(x(i,:),param);
    rem = (1-param)/4;
    w = [rem rem rem rem param];
    cost(i) = Costs(i,:) * w';
end

LCOE = Costs(:,1);
Emissions = Costs(:,2);
DPSP = Costs(:,3);
Dump = Costs(:,4);
REF = ones(dim,1) - Costs(:,5);

% LCOE values on a very different scale so need to plot separately
figure(1)
hold on
plot(params',Emissions,'.-');
plot(params',Dump,'.-');
plot(params',cost,'.-');
plot(params',REF,'.-');
plot(params',LCOE,'.-');
legend('Normalized emissions','Dump','Weighted cost','REF','Normalized LCOE')
xlabel('Value of weight w_5');
ylabel('Objective value [-]');

figure(2)
yyaxis left
plot(params',x(:,3),'.-');
xlabel('Value of weight w_5');
ylabel('BS capacity (kWh)');
yyaxis right
plot(params',DE_hours,'.-');
ylabel('No. of online DE hours per year [h]')
legend('BS capacity', 'DE hours')

%% Oversizing wrt DE rating

RES_sizing_avg = zeros(1,21);
RES_sizing_peak = zeros(1,21);
BS_sizing_daily = zeros(1,21);
params = 0:1:20;

for i = 1:1:21
    RES_rated = x(i,1)*0.255 + x(i,2)*3.5*0.9;
    RES_sizing_avg(i) = RES_rated/8.48;
    RES_sizing_peak(i) = RES_rated/12.52;
    BS_sizing_daily(i) = x(i,3)/203;
end

%plot(params,RES_sizing_avg);
plot(params,BS_sizing_daily);
