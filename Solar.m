function [Ps] = Solar(eta_PV, A_c, I_PV)
% Simple model to calculate ohourly energy output from PV array
% From Tazvinga et al. 2014: Energy dispatch strategy for a photovoltaic–wind–diesel–battery hybrid power system
% A_c = Collector area of the PV array (m^2) 
% eta_PV = Electrical conversion efficiency of the array
% I_PV = Hourly incident solar radiation (kWh/m^2) 
% TODO: Obtain actual I_PV values using GHI data from NASA/NREL (based on lat/lon of chosen location) 

Ps = eta_PV * A_c * I_PV;