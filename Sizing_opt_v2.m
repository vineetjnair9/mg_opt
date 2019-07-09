% Sizing optimization using MATLAB's inbuilt toolboxes

%% Using fminsearch - https://uk.mathworks.com/help/matlab/ref/fminsearch.html
% % But need to constraint variables to be non-negative
% x0 = [1, 1, 10]; % Initial guess for n_s, n_w, Eb_init
% x = fminsearch(@Objective_LI_DE_no_DR,x0)

%% Using https://uk.mathworks.com/help/optim/ug/example-nonlinear-constrained-minimization.html 
% Solve a Constrained Nonlinear Problem, Solver-Based
options = optimoptions(@fmincon,...
    'Display','iter','Algorithm','interior-point');
[x,fval] = fmincon(@Objective_LI_DE_no_DR,[0 0 0],...
    [],[],[],[],[],[],@positive,options)

% Works well to give a reasonable result at least when using 5 equal weights
% But takes a large no. of iterations to converge, might be a local minimum rather than a global one!