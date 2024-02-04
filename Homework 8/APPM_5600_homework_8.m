% Homework 8

%% Problem 3
clc;
clear;
close all;

a = 0; b = 2*pi;
num_x_intervals = 4;
m = 9;
f = @(x) sin(x);
noise = 0*0.05 * randn(m, 1);

[c_2, x_2] = construct_A(a, b, num_x_intervals, m, f, noise);

t_plot_m = 1000;
[t_plot, f_t_plot_2]      = eval_spline(a, b, c_2, x_2, t_plot_m);

figure;
plot(t_plot, f(t_plot), "LineWidth", 2, "DisplayName", "sin(x)");
hold on;
plot(t_plot, f_t_plot_2, "LineWidth", 2, "DisplayName", "Noisy spline");
title("Plot of sin(x) and noisy spline");
xlabel("x");
ylabel("y");
legend;

% Second example
g = @(x) cos(x);

[c_4, x_4] = construct_A(a, b, num_x_intervals, m, g, noise);

[~, g_t_plot_2] = eval_spline(a, b, c_4, x_4, t_plot_m);

figure;
plot(t_plot, g(t_plot), "LineWidth", 2, "DisplayName", "cos(x)");
hold on;
plot(t_plot, g_t_plot_2, "LineWidth", 2, "DisplayName", "Noisy spline");
title("Plot of cos(x) and noisy spline");
xlabel("x");
ylabel("y");
legend;

% Third example
h = @(x) 1./(1+25*x.^(2));

[c_6, x_6] = construct_A(a, b, num_x_intervals, m, h, noise);

[~, h_t_plot_2] = eval_spline(a, b, c_6, x_6, t_plot_m);

figure;
plot(t_plot, h(t_plot), "LineWidth", 2, "DisplayName", "1/(1+25x^{2})");
hold on;
plot(t_plot, h_t_plot_2, "LineWidth", 2, "DisplayName", "Noisy spline");
title("Plot of 1/(1+25x^{2}) and noisy spline");
xlabel("x");
ylabel("y");
legend;