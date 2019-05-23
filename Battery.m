function [soc, Pout] = Battery(soc_prev, P_batt)
% Generalized battery model
% Model obtained from Wu et al. Dynamic economic dispatch of a microgrid: Mathematical models and solution algorithm

delta = 0.01 % Self-discharge rate of battery in %/h
if (P_batt >= 0)
    soc = (1-
    