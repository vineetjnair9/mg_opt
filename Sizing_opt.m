% ERROR: Conversion to double from cvx is not possible
cvx_begin
    variable n_s integer
    variable n_w integer
    variable Eb_init 
    cost = Objective_LI_DE_no_DR_v2([n_s, n_w, Eb_init]);
    minimize cost
    subject to
        x_s >= 1e-4, x_w >= 1e-4, Eb_init >= 1e-4
cvx_end