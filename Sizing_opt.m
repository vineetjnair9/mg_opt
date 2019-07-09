cvx_begin
    variable x_s integer
    variable x_w integer
    variable Eb_init 
    cost = Objective_LI_DE_no_DR(x_s, x_w, Eb_init);
    minimize cost
    subject to
        x_s >= 1e-4, x_w >= 1e-4, Eb_init >= 1e-4
cvx_end