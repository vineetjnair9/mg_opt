%% Plotting out objective function - using Prated
% For various values of inputs -  to show local minima
% Assuming 5 equal weights

Ps_rated_total = 0:1:25; 
Pw_rated_total = 0:1:25;
Eb_init = 0:1:200;
w = [0.2 0.2 0.2 0.2 0.2]; % Weights of different cost terms
x = size(Ps_rated_total,2);
y = size(Pw_rated_total,2);
z = size(Eb_init,2); 

Objectives = zeros(x*y*z,5);
Costs = zeros(x*y*z,1);

tic

for i = 1:1:x
    for j = 1:1:y
        for k = 1:1:z
            for l = 1:1:x*y*z
                input = [Ps_rated_total(i), Pw_rated_total(j), Eb_init(k)];
                Objectives(l,:) = Objective_LI_DE_v6_1h_indPrated(input);
                Costs(l) = Objectives(l,:) * w';
            end
        end
        
    end
end

plotobj_time = toc;

%% Using n_s, n_w