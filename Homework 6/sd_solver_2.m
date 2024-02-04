function [sd_residual_norms, iters, x_old, x_new] = sd_solver_2(A, b, x_0, max_iters, tol)
    iters = 0;
    x_old = x_0;
    sd_residual_norms = [];

    r_k = Inf;
    while (norm(r_k, 2) > tol) && (iters < max_iters)
        r_k = b-A*x_old;
        sd_residual_norms = [sd_residual_norms; norm(r_k, 2)];
        alpha_k = (r_k'*r_k)/(r_k'*A*r_k);
        x_new = x_old + alpha_k*r_k;
        x_old = x_new;
        iters = iters + 1;
    end

    % Compute the next iterate (for problem 2e)
    r_k = b-A*x_old;
    alpha_k = (r_k'*r_k)/(r_k'*A*r_k);
    x_new = x_old + alpha_k*r_k;
end