%% Iterative optimization with feedback 
% Following algorithm of De Paola et al. 
% Can also import custom code as a module into Simulink?

%% (1) Iterative scheme (Flexible demand coordination)
% Initialization phase (1st iteration)
l = 1; % Iteration counter
conv = 0; % Power consumption profiles converged or not?
N = 100; % Number of load agents (individual demands) i.e. price-responsive/flexible appliances
T = 24; % Time interval considered
t = 1:1:T; % Simulation run over 1 day with time discretization step = 1 h (w/ t = 1 being 1st hour starting at 00:00)
P_r = 12e3*ones(1,N); % Assume rated (max) power (in W) = 12 kW for whole population
l_tot = 10000; % Guess for total number of iterations

% Availability window

% Energy requirements for loads are Gaussian normal dist w/ mean = 30 kWh & SD = 1.5 kWh
E_r = 1.5e3*randn(1,N) + 30e3;

% Preallocate array of power profiles for all agents over all times (t) & iteration counts (l) --> Computationally faster
% Or keep concatenating to the multidimensional array?
% But number of iterations needed for convergence not known?

% Feasible power consumption profiles - Enforce by adding extra constraints to the optimization problem
u = ones(N,T,l_tot); % Power profiles for all agents at all times 
% Do as a feasibility problem in CVX - but need to get a set of points and
% not just a single one

% For 1st iteration l=1
for j = 1:N
    u(j,:,1) = 5e3*ones(1,T); % Assume initial power consumption profile to be constant at all times (for all agents)
end

% Assume hourly profile of inflexible demand equal to historical data (0-24h)
% D_i remains fixed at all iterations
% Used historical demand data from the UK as starting example
D_i = [31445;31899;31945;31016;30557;32145;37534;43849;47450;48672;48749;48634;48429;47812;47360;47501;48785;49621;49773;49531;49179;45668;40585;35240]';

% Aggregate demand = sum of flexible and inflexible demands
D = ones(l_tot,T);
D(1,:) = D_i + sum(u(:,1,1));

%% Power scheduling update for (1)
while conv == 0
    conv = 1;
    for j=1:1:N
        l = l+1;
        D(l,:) = D(l-1,:);
        u(j,:,l) = u(j,:,l-1);
        
        for t1=1:1:T
            for t2=1:1:T
                if (t1 ~= t2) && (D(l-1,t1) < D(l-1,t2)) && (u(j,t1,l-1) < P_r(j)) && (u(j,t2,l-1) > 0)
                    % Shift delta of power consumed from times t2 to t1
                    delta = min([P_r(j) - u(j,t1,l-1), u(j,t2,l-1),(D(l-1,t2)-D(l-1,t1))/2]);
                    u(j,t1,l) = u(j,t1,l-1) + delta;
                    D(l,t1) = D(l-1,t1) + delta;
                    u(j,t2,l) = u(j,t2,l-1) - delta;
                    D(l,t2) = D(l-1,t2) + delta;
                    conv = 0;
                end
            end
        end
    end
end

%% One-shot scheme (Flexible demand coordination)
while conv == 0
    conv = 1;
    for j=1:1:N
        l = l+1;
        D(l,:) = D(l-1,:);
        u(j,:,l) = u(j,:,l-1);
        
        for t1=1:1:T
            for t2=1:1:T
                if (t1 ~= t2) && (D(l-1,t1) < D(l-1,t2)) && (u(j,t1,l-1) < P_r(j)) && (u(j,t2,l-1) > 0)
                    % Shift delta of power consumed from times t2 to t1
                    delta = min([P_r(j) - u(j,t1,l-1), u(j,t2,l-1),(D(l-1,t2)-D(l-1,t1))/2]);
                    u(j,t1,l) = u(j,t1,l-1) + delta;
                    D(l,t1) = D(l-1,t1) + delta;
                    u(j,t2,l) = u(j,t2,l-1) - delta;
                    D(l,t2) = D(l-1,t2) + delta;
                    conv = 0;
                end
            end
        end
    end
end

%% Final results
D_nash = D(l,:); % Final aggregate demand profile (converged)
u_nash = u(:,:,l);

figure(1) 
plot(t,D(1,:));
hold on
plot(t,D_nash);
legend('Initial aggregate demand', 'Final aggregate demand')
