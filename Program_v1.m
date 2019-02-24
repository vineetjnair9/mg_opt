%% (1) DESIGN AND PLANNING
%% (a) Sizing storage system
% Function estimates the storage capacity requirements for a specific
% microgrid community
function [out1, out2, out3] = storage(local_demand, ext_power_exhange, local_gen)
    