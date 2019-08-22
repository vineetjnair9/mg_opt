TAC = 17097;
LCOE_grid = 0.125;
CRF = 0.0582;
C_ext = 157470;
load('Load.mat');
P_load = sum(sum(Load));
BEDA_val = (TAC - (LCOE_grid * P_load))/(C_ext*CRF)