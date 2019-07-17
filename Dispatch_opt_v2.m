%% Using fmincon
options = optimoptions(@fmincon,...
    'Display','iter','Algorithm','interior-point');
[x,fval] = fmincon(@Dispatch_obj_LI_DE_no_DR,[zeros(1,24) zeros(1,24)],...
    [],[],[],[],[],[],@dispatch_constraints,options)

%% Using - solve Constrained Nonlinear Optimization, Problem-Based
% https://uk.mathworks.com/help/optim/ug/solve-constrained-nonlinear-optimization-problem-based.html
clear all;
%type Dispatch_obj_LI_DE_no_DR

% Create vector-valued optimization variables
P_DE = optimvar('P_DE', 24)';
P_b_LI = optimvar('P_b_LI', 24)';

% https://uk.mathworks.com/help/optim/ug/fcn2optimexpr.html
[cost, soc_LI, DPSP] = fcn2optimexpr(@Dispatch_obj_LI_DE_no_DR_parallel, P_DE, P_b_LI,'OutputSize',{[1,1],[1,25],[1,1]});
prob = optimproblem;
prob.Objective = cost;

P_DE_rated = 16; % (kW)
SOC_min_LI = 0.1;
SOC_max_LI = 0.9;
P_max_LI = 3.68; % (kW)

Constraint1 = P_DE >= zeros(1,24);
Constraint2 = P_DE <= P_DE_rated*ones(1,24);
Constraint3 = soc_LI >= SOC_min_LI*ones(1,25);
Constraint4 = soc_LI <= SOC_max_LI*ones(1,25);
Constraint5 = P_b_LI <= P_max_LI*ones(1,24);
Constraint6 = P_b_LI >= -P_max_LI*ones(1,24);
Constraint7 = DPSP <= 0.01;

prob.Constraints.Constraint1 = Constraint1;
prob.Constraints.Constraint2 = Constraint2;
prob.Constraints.Constraint3 = Constraint3;
prob.Constraints.Constraint4 = Constraint4;
prob.Constraints.Constraint5 = Constraint5;
prob.Constraints.Constraint6 = Constraint6;
prob.Constraints.Constraint7 = Constraint7;

showproblem(prob)

% Initial guess (start point)
x0.P_DE = ones(1,24);
x0.P_b_LI = ones(1,24);

% But doesn't work since we can't enforce constraints on outputs of
% objective function, only on inputs as optimization variables?

% Trying sqp algorithm, sometimes faster/more accurate than default interior-point
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunEvals',10000);

% Solve the problem
[sol,fval] = solve(prob,x0,'Options',options)