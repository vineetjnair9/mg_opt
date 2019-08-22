%% Using - solve Constrained Nonlinear Optimization, Problem-Based
% https://uk.mathworks.com/help/optim/ug/solve-constrained-nonlinear-optimization-problem-based.html

dim = 24*3;
SOC_initial = linspace(0.1,0.9,9);
Costs_soc = zeros(size(SOC_initial,2),4);
DPSP_soc = zeros(size(SOC_initial,2),1);
P_LI_soc = zeros(size(SOC_initial,2),dim);
P_DE_soc = zeros(size(SOC_initial,2),dim);

% Vary over one week instead of 1 day

for i = 1:1:size(SOC_initial,2)
    % Create vector-valued optimization variables
    P_DE = optimvar('P_DE', dim)';
    P_b_LI = optimvar('P_b_LI', dim)';

    % https://uk.mathworks.com/help/optim/ug/fcn2optimexpr.html
    [Costs, cost, soc_LI, DPSP, P_dump, P_lost, P_load, P_w, P_s, P_RES] = fcn2optimexpr(@Dispatch_obj_LI_DE_v5_varySOC, ...
        P_DE, P_b_LI, SOC_initial(i),'OutputSize',{[1,4],[1,1],[1,dim+1],[1,1],[1,dim],[1,dim],[1,dim],[1,dim],[1,dim],[1,dim]});
    prob = optimproblem;

    prob.Objective = cost;

    P_DE_rated = 8; % (kW)
    SOC_min_LI = 0.1;
    SOC_max_LI = 0.9;
    P_max_LI = 3.68; % (kW)

    % All are bound constraints so can just convert to bounded optimization
    Constraint1 = P_DE >= zeros(1,dim);
    Constraint2 = P_DE <= P_DE_rated*ones(1,dim);
    Constraint3 = soc_LI >= SOC_min_LI*ones(1,dim+1);
    Constraint4 = soc_LI <= SOC_max_LI*ones(1,dim+1);
    Constraint5 = P_b_LI <= P_max_LI*ones(1,dim);
    Constraint6 = P_b_LI >= -P_max_LI*ones(1,dim);
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
    x0.P_DE = 7.*ones(1,dim); % To ensure feasible initial point i.e. DPSP = 0
    x0.P_b_LI = ones(1,dim);

    % Trying sqp algorithm, sometimes faster/more accurate than default interior-point
    options = optimoptions('fmincon','Algorithm','sqp','UseParallel',true,...
        'MaxFunctionEvaluations',100000,'MaxIterations',10000);

    % Solve the problem
    [sol,fval] = solve(prob,x0,'Options',options);
    P_LI_soc(i,:) = sol.P_b_LI;
    P_DE_soc(i,:) = sol.P_DE;
    [Costs_soc(i,:), cost, soc_LI, DPSP_soc(i), P_dump, P_lost, P_load, P_w, P_s, P_RES] = Dispatch_obj_LI_DE_v5_varySOC(P_DE_soc(i,:),P_LI_soc(i,:),SOC_initial(i));
    
end

%%
COE = Costs_soc(:,1);
Emissions = Costs_soc(:,2);
Dump = Costs_soc(:,3);
REF = ones(size(SOC_initial,2),1) - Costs_soc(:,4);
w = [0.2 0.2 0.2 0.2];
cost_soc = Costs_soc * w' + 0.2.*DPSP_soc;

hold on
plot(SOC_initial,COE);
plot(SOC_initial,Emissions);
plot(SOC_initial,DPSP_soc);
plot(SOC_initial,Dump);
plot(SOC_initial,REF);
plot(SOC_initial,cost_soc)
legend('Normalized COE','Emissions','DPSP','Dump energy ratio','REF','Weighted cost')