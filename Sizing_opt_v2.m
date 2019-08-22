 % Sizing optimization using MATLAB's inbuilt toolboxes

%% Using fminsearch - https://uk.mathworks.com/help/matlab/ref/fminsearch.html
% % But need to constrain variables to be non-negative
% x0 = [1, 1, 10]; % Initial guess for n_s, n_w, Eb_init
% x = fminsearch(@Objective_LI_DE_no_DR,x0)

%% Using fmincon (https://uk.mathworks.com/help/optim/ug/example-nonlinear-constrained-minimization.html)
% Solve a Constrained Nonlinear Problem, Solver-Based - fmincon(https://uk.mathworks.com/help/optim/ug/fmincon.html)

% options = optimoptions(@fmincon,...
%     'Display','iter','Algorithm','interior-point','StepTolerance',1e-12,'ConstraintTolerance',1e-12,...
%     'FunctionTolerance',1e-12,'MaxIterations',10000, 'MaxFunctionEvaluations',10000);

options = optimoptions(@fmincon,...
    'Display','iter','Algorithm','interior-point');
[x,fval,exitflag,output] = fmincon(@Objective_LI_DE_v2_1h,[0 0 0],...
    [],[],[],[],[],[],@positive,options)

% Works well to give a reasonable result at least when using 5 equal weights
% But takes a large no. of iterations to converge, might be a local minimum rather than a global one!

%% GlobalSearch (https://uk.mathworks.com/help/gads/globalsearch.html)
% No longer need to specify constraints on inputs explicitly since it's now
% a bounded constrained optimization (positive constraint is enforced by lower bound)
% Maybe remove upper bound constraint - don't need to explicitly specify this?

rng default % For reproducibility
ms = MultiStart('UseParallel',true);
gs = GlobalSearch(ms);
options = optimoptions(@fmincon, 'Display','iter','Algorithm','active-set',...
    'MaxIterations',2000,'MaxFunctionEvaluations',2000,'UseParallel',true);
    %'PlotFcn',{@gsplotbestf,@gsplotfunccount});
% Can't plot with GS or MS while using parallel processing
problem = createOptimProblem('fmincon','objective',@Objective_LI_DE_v2_1h,'x0',[50 3.5 50],...
    'lb',[0,0,0],'ub',[100,30,200],'options',options);
tic
[x, fval] = run(gs,problem)
gs_time = toc; % Time taken for optimization

%% Multiple starting point search (https://uk.mathworks.com/help/gads/multistart.html)
rng default % For reproducibility
options = optimoptions(@fmincon,'Display','iter',...
    'Algorithm','active-set','MaxIterations',2000,'MaxFunctionEvaluations',2000);
problem = createOptimProblem('fmincon','objective',...
    @Objective_LI_DE_v2_1h,'x0',[50 3.5 50],...
    'lb',[0,0,0],'ub',[100,30,200],'options',options);
ms = MultiStart('UseParallel',true);
tic
[x,f] = run(ms,problem,85)
ms_time = toc;

%% Direct pattern search (https://uk.mathworks.com/help/gads/patternsearch.html)
% Pattern search using several different initial/start points 
rng default % For reproducibility
lb = [0,0,0];
ub = [100,30,200];
A = [];
b = [];
Aeq = [];
beq = [];

x = zeros(85,3);
fval = zeros(85,1);
exitflag = zeros(85,1);
options = optimoptions('patternsearch','UseParallel',true,'UseCompletePoll',true,...
    'UseVectorized',false);

%PlotFcn',{@psplotbestf,psplotbestx,psplotmeshsize,psplotfuncount}

tic
for i = 1:1:85
    x0 = lb + rand(size(lb)).*(ub - lb);
    [x(i,:), fval(i), exitflag(i)] = patternsearch(@Objective_LI_DE_v2_1h,x0,A,b,Aeq,beq,lb,ub,options);
end

[min_obj, index] = min(fval)
x_opt = x(index,:)
ps_time = toc;

%% Genetic algorithm (GA) 
rng default % For reproducibility
lb = [0,0,0];
ub = [100,30,200];
fitness_function = @Objective_LI_DE_v2_1h;
nvars = 3; % number of variables
options = optimoptions('ga','UseParallel', true, 'UseVectorized', false,...
    'PlotFcn','gaplotgenealogy');
%@gaplotbestf,@gaplotbestindiv,@gaplotrange,
tic
[x,fval] = ga(fitness_function,nvars,[],[],[],[],lb,ub,[],options)
ga_time = toc;

%% Plotting objective
plotobjective(@Objective_LI_DE_v2_1h,[0 30; 0 10; 0 100]);

%% Genetic algorithm (GA) - vectorized
rng default % For reproducibility
lb = [0,0,0];
ub = [100,20,200];
fitness_function = @Objective_LI_DE_v2_1h_vectorized;
nvars = 3; % number of variables
options = optimoptions('ga','UseParallel', true, 'UseVectorized', true);
tic
[x,fval] = ga(fitness_function,nvars,[],[],[],[],lb,ub,[],options)
ga_time = toc;

% But not easy to implement in this case due to the manner in which the 
% objective function is written, preallocate arrays for all the powers etc.
%% PSO
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false,'PlotFcn','pswplotbestf');
nvars = 3;
lb = [0,0,0];
ub = [100,30,200];
tic
[x,fval,exitflag] = particleswarm(@Objective_LI_DE_v2_1h,nvars,lb,ub,options)
pso_time = toc;

%% Surrogate optimization
rng default
lb = [0,0,0];
ub = [100,30,200];
options = optimoptions('surrogateopt','UseParallel',true,'MaxFunctionEvaluations',...
    2000,'PlotFcn','surrogateoptplot');
tic
[x,fval] = surrogateopt(@Objective_LI_DE_v2_1h,lb,ub,options);
sur_time = toc;

%% Simulated annealing
rng default;
lb = [0,0,0];
ub = [100,30,200];

x0 = [50 3.5 50]; % initial point
options = optimoptions('simulannealbnd','PlotFcns',...
    {@saplotbestf,@saplotf,@saplotstopping});
tic
[x,~,exitFlag,output] = simulannealbnd(@Objective_LI_DE_v2_1h,x0,lb,ub,options)
sa_time = toc;

%% PSO - but optimizing wrt rated powers instead
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false,'PlotFcn','pswplotbestf');
nvars = 3;
lb = [0,0,0];
ub = [25,25,200];
tic
[x,fval,exitflag] = particleswarm(@Objective_LI_DE_v8_1h_optPrated,nvars,lb,ub,options)
pso_time = toc;

%% Pareto set using gamultiobj
rng default;
lb = [0,0,0];
ub = [100,10,200]; % for n_s, n_w, Eb_init
tic
options = optimoptions('gamultiobj','UseParallel',true,'UseVectorized',false,...
    'ParetoFraction',0.50);
[x_p,fval] = gamultiobj(@sizing_pareto,3,[],[],[],[],lb,ub,options);
gamultiobj_time = toc;

%% Pareto set with rated total powers as inputs instead of n_s/n_w
rng default;
options = optimoptions('gamultiobj','UseParallel',true,...
    'UseVectorized',false,'ParetoFraction',0.50);
lb = [0,0,0];
ub = [25,25,200];
tic
[x_p,fval] = gamultiobj(@sizing_pareto,3,[],[],[],[],lb,ub,[],options);
gamultiobj_time = toc;

%% Pareto set/front using paretosearch (pattern search algorithm)
rng default;
options = optimoptions('paretosearch','UseParallel',true,...
    'UseVectorized',false,'ParetoSetSize',200);
lb = [0,0,0];
ub = [100,10,200]; % for rated P_s, P_w, Eb_init
tic
[x_p,fval] = paretosearch(@sizing_pareto,3,[],[],[],[],lb,ub,[],options);
paretosearch_time = toc;

%% Pareto set with rated total powers as inputs instead of n_s/n_w
rng default;
options = optimoptions('paretosearch','UseParallel',true,...
    'UseVectorized',false,'ParetoSetSize',200);
lb = [0,0,0];
ub = [25,25,200];
tic
[x_p,fval] = paretosearch(@sizing_pareto,3,[],[],[],[],lb,ub,[],options);
paretosearch_time = toc;

%%
w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
cost = fval * w';

%% 3-D pareto set plot in control variable (inputs) space
% 3D scatter plot 
% Plotting optimal values of n_s, n_w and Eb_init
scatter3(x_p(:,1),x_p(:,2),x_p(:,3),10,'b','filled')
xlabel('No. of solar PV arrays');
ylabel('No. of wind turbines');
zlabel('Battery capacity (kWh)');

%% Interpolated surface plot
figure(1)
x = x_p(:,1);
y = x_p(:,2);
z = x_p(:,3);
F = scatteredInterpolant(x,y,z,'natural','none');
sgr = linspace(min(x),max(x));
ygr = linspace(min(y),max(y));
[XX,YY] = meshgrid(sgr,ygr);
ZZ = F(XX,YY);
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,5,'k','filled')
%scatter3(x,y,z,'k.')

%% Griddata
figure(2)
x = x_p(:,1);
y = x_p(:,2);
z = x_p(:,3);
sgr = linspace(min(x),max(x));
ygr = linspace(min(y),max(y));
[XX,YY] = meshgrid(sgr,ygr);
ZZ = griddata(x,y,z,XX,YY,'natural');
%mesh(xq,yq,zq);
surf(XX,YY,ZZ,'LineStyle','none');
hold on
plot3(x,y,z,'k.')

%% Rotated views
figure
subplot(2,2,1)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
subplot(2,2,2)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-148,8)
subplot(2,2,3)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-180,8)
subplot(2,2,4)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-300,8)

%% Rotated separate figures
font = 15;
figure(1)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
xlabel('No. of solar PV arrays');
ylabel('No. of wind turbines');
zlabel('Battery capacity (kWh)');
set(gca,'FontSize',font);
figure(2)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-148,8)
xlabel('No. of solar PV arrays');
ylabel('No. of wind turbines');
zlabel('Battery capacity (kWh)');
set(gca,'FontSize',font);
figure(3)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-180,8)
xlabel('No. of solar PV arrays');
ylabel('No. of wind turbines');
zlabel('Battery capacity (kWh)');
figure(4)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-300,8)
xlabel('No. of solar PV arrays');
ylabel('No. of wind turbines');
zlabel('Battery capacity (kWh)');
set(gca,'FontSize',font);

%% 
xlabel('PV capacity (kW)');
ylabel('WT capacity (kW)');
zlabel('BS capacity (kWh)');

%% Rotated separate figures
figure(1)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
xlabel('PV capacity (kW)');
ylabel('WT capacity (kW)');
zlabel('BS capacity (kWh)');
figure(2)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-148,8)
xlabel('PV capacity (kW)');
ylabel('WT capacity (kW)');
zlabel('BS capacity (kWh)');
figure(3)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-180,8)
xlabel('PV capacity (kW)');
ylabel('WT capacity (kW)');
zlabel('BS capacity (kWh)');
figure(4)
surf(XX,YY,ZZ,'LineStyle','none')
hold on
scatter3(x,y,z,'k.');
hold off
view(-300,8)
xlabel('PV capacity (kW)');
ylabel('WT capacity (kW)');
zlabel('BS capacity (kWh)');

%% 2D Pareto front curves for inputs
n_s = x_p(:,1);
n_w = x_p(:,2);
Eb_init = x_p(:,3);

figure(1)
scatter(n_s,n_w,10,'b','filled');
xlabel('PV capacity (kW)');
ylabel('WT capacity (kW)');

figure(2)
scatter(n_s,Eb_init,10,'b','filled');
xlabel('PV capacity (kW)');
ylabel('Battery capacity (kWh)');

figure(3)
scatter(n_w,Eb_init,10,'b','filled');
xlabel('WT capacity (kW)');
ylabel('Battery capacity (kWh)');

%% Pareto set for objectives
LCOE = fval(:,1); 
Emissions = fval(:,2); 
DPSP = fval(:,3); 
Dump = fval(:,4);
%REF = ones(size(fval,1),1) - fval(:,5);

scatter(Dump,Emissions,5,'b','filled')
xlabel('Dump energy ratio');
ylabel('Normalized emissions');

%% fsolve to obtain pareto front of inputs

Cost = @ (x) Objective_LI_DE_v2_1h(x) - 0.1367;
OPTIONS = optimoptions('fsolve','Algorithm','Levenberg-Marquardt',...
    'StepTolerance',1e-10);
[X] = fsolve(Cost,[10 3 100],OPTIONS);

















