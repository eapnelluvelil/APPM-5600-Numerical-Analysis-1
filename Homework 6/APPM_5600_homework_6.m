%% Homework 6

%% Problem 2
clear;
clc;
close all;

fprintf("Problem 2\n");
fprintf("---------\n\n");

rng(0);

max_iters = 1000;
tol = 1e-10;
n = 500;
taus = [0.01, 0.05, 0.1, 0.2];
x_0 = zeros(n, 1);
b = rand(n, 1);

% Create a matrix with entries from the uniform distribution on [-1, 1]
% and extract the strictly upper triangular portion of it
U_strict = triu(-1 + (1-(-1))*rand(n), 1);

diag_A = diag(ones(n, 1));

f1 = figure;
f2 = figure;

for i=1:length(taus)
    U_tau = U_strict;
    U_tau(abs(U_tau) > taus(i)) = 0;
    A_tau = diag_A + U_tau + U_tau';

    % Compute the true solution using backslash
    x_star = A_tau\b;

    fprintf("Smallest eigenvalue of A_tau: %e\n\n", min(eig(A_tau)));

    % [sd_residual_norms, sd_iters] = sd_solver(A_tau, b, x_0, max_iters, tol);
    [sd_residual_norms, sd_iters, x_old_sd, x_new_sd] = sd_solver_2(A_tau, b, x_0, max_iters, tol);

    fprintf("tau = %0.2f, norm of last residual: %e, iterations of SD: %d\n", ...
        taus(i), sd_residual_norms(end), sd_iters);

    % Compute error bound for steepest descent iteration
    e_old_sd = norm((x_old_sd - x_star)'*A_tau*(x_old_sd - x_star));
    e_new_sd = norm((x_new_sd - x_star)'*A_tau*(x_new_sd - x_star));
    A_tau_eigs = eig(A_tau);
    sd_error_const = (max(A_tau_eigs)-min(A_tau_eigs))/(max(A_tau_eigs)+min(A_tau_eigs));
    fprintf("Error bound for steepest descent at next iterate: %e <= %e\n\n", ...
            e_new_sd, sd_error_const*e_old_sd);
    
    figure(f1);
    semilogy(1:sd_iters, sd_residual_norms, "LineWidth", 2, ...
             "DisplayName", strcat("tau = ", string(taus(i))));
    hold on;
    
    [cg_residual_norms, cg_iters, x_old_cg] = cg_solver(A_tau, b, x_0, max_iters, tol);
    
    fprintf("tau = %0.2f, norm of last residual: %e, iterations of CG: %d\n", ...
        taus(i), cg_residual_norms(end), cg_iters);

    % Compute error bound for conjugate gradient iteration
    cond_A_tau = cond(A_tau);
    cg_error_const = 2*((1-sqrt(1/cond_A_tau))/(1+sqrt(1/cond_A_tau)))^(cg_iters);
    e_new_cg = norm((x_old_cg - x_star)'*A_tau*(x_old_cg - x_star));
    x_star_A_norm = norm(x_star'*A_tau*x_star);
    fprintf("Error bound for conjugate gradient at latest iterate: %e <= %0.16e\n\n", ...
            e_new_cg, cg_error_const*x_star_A_norm);
    
    figure(f2);
    semilogy(1:cg_iters, cg_residual_norms, "LineWidth", 2, ...
             "DisplayName", strcat("tau = ", string(taus(i))));
    hold on;
end

figure(f1);
title("Number of iterations versus 2-norm of residuals for steepest descent iterations for various tau");
xlabel("n");
ylabel("||r_{n}||_{2}");
legend;

figure(f2);
title("Number of iterations versus 2-norm of residuals for conjugate gradient iterations for various tau");
xlabel("n");
ylabel("||r_{n}||_{2}");
legend;

%% Problem 4
clear;
close all;

fprintf("\n\nProblem 4\n");
fprintf("---------\n\n")

tol = 1e-7;

g = @(x, y) [x; y] - [0.016 -0.17; 0.52 -0.26]*[3*x.^(2) + 4*y.^(2) - 1; y.^(3) - 8*x.^(3) - 1];

x_old = [-0.5, 0.25];
x_new = g(x_old(1), x_old(2));
fixed_pt_iters = 1;

while norm(x_new-x_old, 2)/norm(x_old, 2) > tol
    x_old = x_new;
    x_new = g(x_old(1), x_old(2));
    fixed_pt_iters = fixed_pt_iters + 1;
end

fprintf("Converged fixed point iterate: (%f, %f)\n", x_new(1), x_new(2));
fprintf("Number of fixed point iterations to get relative error below %e: %d\n", ...
        tol, fixed_pt_iters);