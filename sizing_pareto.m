function f = sizing_pareto(x)

Costs = Objective_LI_DE_v6_1h_indPrated(x);
f(1) = Costs(1);
f(2) = Costs(2);
f(3) = Costs(3);
%f(4) = Costs(4);
%f(5) = Costs(5);

end