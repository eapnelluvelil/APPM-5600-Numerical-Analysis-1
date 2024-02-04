function [sd_residual_norms, iters] = sd_solver(A, b, x_0, max_iters, tol)
iters = 1;
r_k = b-A*x_0;
sd_residual_norms = [norm(r_k, 2)];
alpha_k = (r_k'*r_k)/(r_k'*A*r_k);
x_k = x_0 + alpha_k*r_k;

while (norm(r_k, 2) > tol) && (iters < max_iters)
    r_k = b-A*x_k;
    sd_residual_norms = [sd_residual_norms; norm(r_k, 2)];
    alpha_k = (r_k'*r_k)/(r_k'*A*r_k);
    x_k = x_k + alpha_k*r_k;
    iters = iters + 1;
end
end