%% LI
Eb_init = 100;
E_C = zeros(1,5475);
E_C(1) = Eb_init; 
cycles_LI = 1;
cycles = 0:1:5474;

for i = 2:1:5475
    E_C(i) = Eb_init*(1 - (cycles_LI*0.055e-3));
    cycles_LI = cycles_LI+1;
end

plot(cycles,E_C);