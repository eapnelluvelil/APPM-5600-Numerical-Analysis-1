%% Homework 7

%% Problem 2
clc;
clear;
close all;

f = @(x, y) exp(x).*sin(y);

x_interp = [0; 0; 1; 1; 2; 2];
y_interp = [0; 2; 0; 2; 1; 3];

A = [ones(size(x_interp)), x_interp, y_interp, x_interp.*y_interp, ...
     x_interp.^(2), y_interp.^(2)];

b = f(x_interp, y_interp);

% Solve for the multivariate interpolant coefficients
c = A\b;
fprintf("Coefficients for question 2b: \n");
disp(c);

f_interp = @(c, x, y) c(1)*ones(size(x)) + c(2)*x + c(3)*y + c(4)*x.*y + c(5)*x.^(2) + c(6)*y.^(2);

a = -1; b = 3; n = 25;
[XX, YY] = meshgrid(linspace(a, b, n), linspace(a, b, n));
ZZ_f = f(XX, YY);
ZZ_f_interp = f_interp(c, XX, YY);

figure(1);
surf_f = surf(XX, YY, ZZ_f);
xlabel("x");
ylabel("y");
zlabel("z");
title("Surface plot of f(x) = e^{x} sin(y)");
legend(surf_f, {'f'});

figure(2);
surf_f_interp = surf(XX, YY, ZZ_f_interp);
xlabel("x");
ylabel("y");
zlabel("z");
title("Surface plot of p_{6} interpolant of f(x) = e^{x} sin (y)");
legend(surf_f_interp, {'p_6'});