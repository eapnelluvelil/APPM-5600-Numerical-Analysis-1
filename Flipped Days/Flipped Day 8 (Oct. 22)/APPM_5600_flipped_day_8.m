%% Flipped Day 8 (October 22, 2021)
clc;
clear;
close all;

% Problem 1
f = @(x) exp(sin(x));

n = 10;
x = (0:1:(2*n))'*((2*pi)/(2*n+1));
z = exp(1i*x);

V = zeros(2*n+1);
for j = 1:(2*n+1)
    V(:, j) = z.^(j-1);
end


f_bar = (z.^(n)).*(f(x));
c = (V'*f_bar)/(2*n+1);

n_fine = n;
x_fine = (0:1:(2*n_fine))'*((2*pi)/(2*n_fine+1));
z_fine = exp(1i*x_fine);

figure;
plot(x_fine, f(x_fine));
hold on;
plot(x_fine, eval_trig_interp(c, z_fine));


function val = eval_trig_interp(c, z)
    V = zeros(length(z), length(c));
    for j = 1:length(c)
        V(:, j) = z.^(j);
    end
    val = V*c;
end