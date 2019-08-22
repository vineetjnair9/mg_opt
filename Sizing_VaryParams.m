%% Varying parameters in sizing optimization
% Including variable weights on multiobjective
% Optimizing w.r.t n_s, n_w, Eb_init

%% Varying DG rated power
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false);
nvars = 3;
lb = [0,0,0];
ub = [100,30,200];

% Use opt. x solutions later to calculate 5 individual objectives, DE hours & BS cycles!

DE_rating = 0:1:20; % in kW
x = zeros(size(DE_rating,2),3);
fval = zeros(size(DE_rating,2),1);

tic
for i = 1:1:size(DE_rating,2)
    P_DE_rated = DE_rating(i);
    Objective = @(x) Objective_LI_DE_v9_1h_ind_VaryParams(x,P_DE_rated);
    [x(i,:),fval(i)] = particleswarm(Objective,nvars,lb,ub,options);
end
param_time = toc;

% Normalized LCOE and emissions values for lower ratings may be artificially inflated
% slightly, since these are now being considered relative to a smaller size
% DG that may not satisfy all load in base case

%% Varying fuel (diesel price)
% Can also account for effect of carbon tax or fossil fuel subsidies
% Range of global diesel prices - https://www.globalpetrolprices.com/diesel_prices/
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false);
nvars = 3;
lb = [0,0,0];
ub = [100,30,200];

C_DE = linspace(0.2,12,20); % in $/gal
x = zeros(size(C_DE,2),3);
fval = zeros(size(C_DE,2),1);

tic
for i = 1:1:size(C_DE,2)
    DE_price = C_DE(i);
    Objective = @(x) Objective_LI_DE_v9_1h_ind_VaryParams(x,DE_price);
    [x(i,:),fval(i)] = particleswarm(Objective,nvars,lb,ub,options);
end
param_time = toc;

%% Varying nominal interest/discount rate
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false);
nvars = 3;
lb = [0,0,0];
ub = [100,30,200];

inflation = linspace(0,0.2,20); % in (%)
x = zeros(size(inflation,2),3);
fval = zeros(size(inflation,2),1);

tic
for i = 1:1:size(inflation,2)
    interest_rate = inflation(i);
    Objective = @(x) Objective_LI_DE_v9_1h_ind_VaryParams(x,interest_rate);
    [x(i,:),fval(i)] = particleswarm(Objective,nvars,lb,ub,options);
end
param_time = toc;

%% Varying inflation rate
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false);
nvars = 3;
lb = [0,0,0];
ub = [100,30,200];

inflation = linspace(0,0.15,20); % in (%)
x = zeros(size(inflation,2),3);
fval = zeros(size(inflation,2),1);

tic
for i = 1:1:size(inflation,2)
    inflation_rate = inflation(i);
    Objective = @(x) Objective_LI_DE_v9_1h_ind_VaryParams(x,inflation_rate);
    [x(i,:),fval(i)] = particleswarm(Objective,nvars,lb,ub,options);
end
param_time = toc;

%% Vary LI BS price
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false);
nvars = 3;
lb = [0,0,0];
ub = [100,30,200];

BS_price = linspace(50,300,20); % in ($/kWh)
x = zeros(size(BS_price,2),3);
fval = zeros(size(BS_price,2),1);

tic
for i = 1:1:size(BS_price,2)
    BS = BS_price(i);
    Objective = @(x) Objective_LI_DE_v9_1h_ind_VaryParams(x,BS);
    [x(i,:),fval(i)] = particleswarm(Objective,nvars,lb,ub,options);
end
param_time = toc;

%% Varying relative values of weights
rng default % For reproducibility of results (needed for random algorithm)
options = optimoptions('particleswarm', 'UseParallel', true, ...
    'UseVectorized', false);
nvars = 3;
lb = [0,0,0];
ub = [100,30,200];

w3 = linspace(0,1,11); 
dim = size(w3,2);
x = zeros(dim,3);
fval = zeros(dim,1);

tic
for i = 1:1:dim
    w = w3(i);
    Objective = @(x) Objective_LI_DE_v9_1h_ind_VaryParams(x,w);
    [x(i,:),fval(i)] = particleswarm(Objective,nvars,lb,ub,options);
end
param_time = toc;









