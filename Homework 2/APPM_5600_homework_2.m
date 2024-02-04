%% Problem 2
fprintf("Homework 2, Problem 2\n");
% To call our bisection implementation
addpath(genpath("C:\Users\eappe\Google Drive\1st Year\APPM 5600\Homework 1"));
close all;

T_i = 20;
T_s = -15;
alpha = 0.138 * 1e-6;
tol = 1e-13;

% x is in meters, t is in seconds
% 60 days after the cold snap, which is 60*24*60*60 seconds
t = 60*24*60*60; 
% Temperature function
T = @(x) (T_i - T_s)*erf(x./(2*sqrt(alpha*t))) + T_s; 
% Derivative of temperature function
T_prime = @(x) (T_i - T_s)*(1/sqrt(pi*alpha*t))*exp(-((x.^(2))./4*alpha*t)); 

a =  0;
b =  1;
x = linspace(a, b, 1000);

figure;
plot(x, T(x), "LineWidth", 2);
xlabel("x (meters below surface)");
ylabel("T(x, t) (temperature (in degrees Celsius))");
title("Plot of T(x, t)");

% Find a root of T on the interval [x_0, x_bar] via bisection
[r, iters] = bisect(a, b, T, tol);
fprintf("Root computed via bisection: %0.16f (%d iterations)\n", r, iters);

% Find a root of T on the interval [x_0, x_bar] via Newton's method
% Start off with an initial guess close to 0
x_0 = 0.01;
n_max = 100;
[r, iters, ~] = newton(x_0, T, T_prime, tol, n_max); 
fprintf("Root computed via Newton's method: %0.16f (%d iterations, x_0 = %0.16f)\n", ...
        r(iters), iters, x_0);

% Use Newton's method again, but make the initial guess the endpoint of the
% interval [a, b], i.e., x_0 = b
x_0 = b;
[r, iters, ~] = newton(x_0, T, T_prime, tol, n_max); 
fprintf("Root computed via Newton's method: %0.16f (%d iterations, x_0 = %0.16f)\n", ...
        r(iters), iters, x_0);

%% Problem 4
fprintf("Homework 2, Problem 4\n");
close all;

p = 5;
f = @(x) ((x-1).^(p)).*exp(x);
f_prime = @(x) (p*(x-1).^(p-1)).*exp(x) + ((x-1).^(p)).*exp(x);
% f_prime = @(x) ((x-1).^(4)).*(x+4).*exp(x);

tol = 1e-15;
n_max = 150;

x_0 = 0;
[r, iters, ~] = newton(x_0, f, f_prime, tol, n_max);
fprintf("Root computed via Newton's method: %0.16f (%d iterations, x_0 = %0.16f)\n", ...
        r(iters), iters, x_0);
fprintf("Relative error between approximate root and 1: %0.16e\n\n", ...
        abs(r(iters) - 1));

Newton_iterates = r(1:iters);
Newton_relative_errors = abs(r(1:iters) - 1);
T1 = table(Newton_iterates, Newton_relative_errors);
disp(T1);

figure;
semilogy(1:iters, abs(r(1:iters) - 1), "LineWidth", 2);
xlabel("k");
ylabel("log(|e_{k}|)");
title("Plot of relative errors computed via Newton's method");

[r, iters, ~] = modified_newton(x_0, f, f_prime, p, tol, n_max);
fprintf("Root computed via modified Newton's method: %0.16f (%d iterations, x_0 = %0.16f)\n", ...
        r(iters), iters, x_0);
fprintf("Relative error between approximate root and 1: %0.16e\n\n", ...
    abs(r(iters) - 1));

Modified_Newton_iterates = r(1:iters);
Modified_Newton_relative_errors = abs(r(1:iters) - 1);
T2 = table(Modified_Newton_iterates, Modified_Newton_relative_errors);
disp(T2);

figure;
semilogy(1:iters, abs(r(1:iters) - 1), "LineWidth", 2);
xlabel("k");
ylabel("log(|e_{k}|)");
title("Plot of relative errors computed via modified Newton's method");