function [cg_residual_norms, iters, x_old] = cg_solver(A, b, x_0, max_iters, tol)
iters = 1;
x_old = x_0;
r_old = b-A*x_old;
cg_residual_norms = [norm(r_old, 2)];
p_old = r_old;

while (norm(r_old, 2) > tol) && (iters < max_iters)
    alpha_new = (r_old'*p_old)/(p_old'*A*p_old);
    x_new = x_old + alpha_new*p_old;
    r_new = r_old - alpha_new*A*p_old;
    cg_residual_norms = [cg_residual_norms; norm(r_new, 2)];
    beta_new = (r_new'*r_new)/(r_old'*r_old);
    p_new = r_new + beta_new*p_old;
    iters = iters + 1;
    r_old = r_new;
    p_old = p_new;
    x_old = x_new; 
end
end