%% Problem 2
fprintf("Homework 1, Problem 2\n");
close all;

x = 1.920:0.001:2.080;

p_direct = (x-2).^9;
p_coeffs = [1, -18, 144, -672, 2016, -4032, 5376, -4608, 2304, -512];
p_coeffs = polyval(p_coeffs, x);

% Plot the polynomial p, which is evaluated directly and via the expanded
% form
figure;
plot(x, p_direct, "LineWidth", 2);
hold on;
plot(x, p_coeffs, "LineWidth", 2);
legend("p evaluated directly", "p evaluated via coefficients");
xlabel("x");
ylabel("p");
title("Plot of p_{direct} and p_{expanded} for x near 2");

%% Problem 3
fprintf("Homework 1, Problem 3\n");
close all;

delta = logspace(-16, 0);

% Various evaluations of cos(x + delta) - cos(x)
cos_diff_direct = @(x, delta) cos(x + delta) - cos(x);
cos_diff_identity = @(x, delta) -2*sin(x + delta/2).*sin(delta/2);
cos_diff_taylor = @(x, delta) -delta.*sin(x) - ...
    (1/2)*((delta).^(2)).*cos(x+delta/2);

% Compute the relative errors between the approximations to 
% cos(x + delta) - cos(x)
rel_err_direct_identity = @(x) abs((cos_diff_direct(x, delta) - ...
    cos_diff_identity(x, delta))./cos_diff_direct(x, delta));
rel_err_direct_taylor = @(x) abs((cos_diff_direct(x, delta) - ...
    cos_diff_taylor(x, delta))./cos_diff_direct(x, delta));

% Compare the two ways of evaluating cos(x+delta) - cos(x) when x is large
% and small
x1 = 1e10;
x2 = 4*atan(1);

figure;
loglog(delta, rel_err_direct_identity(x1), "LineWidth", 2);
hold on;
loglog(delta, rel_err_direct_taylor(x1), "LineWidth", 2);
hold on;
loglog(delta, rel_err_direct_identity(x2), "LineWidth", 2);
hold on;
loglog(delta, rel_err_direct_taylor(x2), "LineWidth", 2);
legend(compose("Trig identity (for x = %1.2e)", x1), compose("Taylor expansion (for x = %1.2e)", x1), ...
       compose("Trig identity (for x = %f)", x2), compose("Taylor expansion (for x = %f)", x2));
xlabel("$\delta$", "Interpreter", "latex");
ylabel("$\log$ (Relative error)", "Interpreter", "latex");
title("Relative error when $\cos (x + \delta) - \cos (x)$ is evaluated two ways", ...
      "Interpreter", "latex");

%% Problem 6
fprintf("Homework 1, Problem 6\n");
close all;

a = 4.8;
b = 5.31;
tol = 1e-4;

f_direct = @(x) (x-5).^(9);

fc = [1; -45; 900; -10500; 78750; -393750; 1312500; -2812500; 3515625; -1953125];
f_ceoffs_eval = @(x) polyval(fc, x);

[r, iters_1] = bisect(a, b, f_direct, tol);
fprintf("root = %0.16f, f(r) = %0.16f, number of iterations: %d\n", r, f_direct(r), iters_1);

[r, iters_2] = bisect(a, b, f_coeffs_eval, tol);
fprintf("root = %0.16f, f(r) = %0.16f, number of iterations: %d\n", r, f_coeffs_eval(r), iters_2);

figure;
x = [a:0.001:b];
plot(x, abs(f_direct(x)-f_coeffs_eval(x)), "LineWidth", 2);
xlabel("x");
ylabel("|f_{direct}(x) - f_{expanded}(x)|");
title("Plot of |f_{direct} - f_{expanded}| for x near 5");
