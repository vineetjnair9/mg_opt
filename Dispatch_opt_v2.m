%% Using - solve Constrained Nonlinear Optimization, Problem-Based
% https://uk.mathworks.com/help/optim/ug/solve-constrained-nonlinear-optimization-problem-based.html

% Create vector-valued optimization variables
P_DE = optimvar('P_DE', 24)';
P_b_LI = optimvar('P_b_LI', 24)';

% https://uk.mathworks.com/help/optim/ug/fcn2optimexpr.html
[cost, Costs, soc_LI, DPSP, P_dump, P_lost, P_load, P_w, P_s, P_RES] = fcn2optimexpr(@Dispatch_obj_LI_DE_v3, P_DE, P_b_LI,'OutputSize',{[1,1],[1,25],[1,1],...
    [1,24],[1,24],[1,24],[1,24],[1,24],[1,24]});
prob = optimproblem;

prob.Objective = cost;

P_DE_rated = 8; % (kW)
SOC_min_LI = 0.1;
SOC_max_LI = 0.9;
P_max_LI = 3.68; % (kW)

% All are bound constraints so can just convert to bounded optimization
Constraint1 = P_DE >= zeros(1,24);
Constraint2 = P_DE <= P_DE_rated*ones(1,24);
Constraint3 = soc_LI >= SOC_min_LI*ones(1,25);
Constraint4 = soc_LI <= SOC_max_LI*ones(1,25);
Constraint5 = P_b_LI <= P_max_LI*ones(1,24);
Constraint6 = P_b_LI >= -P_max_LI*ones(1,24);
Constraint7 = DPSP <= 0.0100;

prob.Constraints.Constraint1 = Constraint1;
prob.Constraints.Constraint2 = Constraint2;
prob.Constraints.Constraint3 = Constraint3;
prob.Constraints.Constraint4 = Constraint4;
prob.Constraints.Constraint5 = Constraint5;
prob.Constraints.Constraint6 = Constraint6;
prob.Constraints.Constraint7 = Constraint7;

showproblem(prob)

% Initial guess (start point)
x0.P_DE = 7.*ones(1,24); % To ensure feasible initial point i.e. DPSP = 0
x0.P_b_LI = ones(1,24);

% Trying sqp algorithm, sometimes faster/more accurate than default interior-point
options = optimoptions('fmincon','Algorithm','sqp','UseParallel',true,...
    'MaxFunctionEvaluations',100000,'MaxIterations',10000,'Display','iter');

% Solve the problem
tic
[sol,fval] = solve(prob,x0,'Options',options)
dispatch_time = toc;