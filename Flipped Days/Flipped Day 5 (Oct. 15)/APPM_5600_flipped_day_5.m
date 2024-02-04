clc;
clear;
close all;

f = @(x) x./(1+25*x.^(2));
f_prime = @(x) (1-25*x.^(2))./((1+25*x.^(2)).^(2));

a = -1;
b = 1;

n_grid = 1000;
grid_pts = linspace(a, b, n_grid);

% Construct Lagrange interpolant using equispaced points
n_equi = 21;
equi_pts_lagrange = linspace(a, b, n_equi);
% equi_pts = cos(((2*[1:n_equi]-1)/(2*n_equi))*pi);
f_equi_pts = f(equi_pts_lagrange);

lagrange_interp = zeros(size(grid_pts));
for i = 1:length(equi_pts_lagrange)
    l_i = ones(size(grid_pts));
    for j = 1:length(equi_pts_lagrange)
        if j ~= i
            l_i = l_i .* ((grid_pts-equi_pts_lagrange(j))/(equi_pts_lagrange(i)-equi_pts_lagrange(j)));
        end
    end
    f_l_i = l_i*f_equi_pts(i);
    lagrange_interp = lagrange_interp + f_l_i;
end

figure(1);
plot(grid_pts, f(grid_pts));
hold on;
plot(grid_pts, lagrange_interp);

% Construct Hermite interpolant using equispaced points
n_equi_hermite = 11;
equi_pts_hermite = linspace(a, b, n_equi_hermite);
f_equi_pts_hermite = f(equi_pts_hermite);
f_prime_equi_pts_hemrite = f_prime(equi_pts_hermite);

hermite_interp = zeros(size(grid_pts));
for i = 1:length(equi_pts_hermite)
    l_i = ones(size(grid_pts));
    h_i = ones(size(grid_pts));

    for j = 1:length(equi_pts_hermite)
        if j ~= i
            l_i = l_i .* ((grid_pts-equi_pts_hermite(j))/(equi_pts_hermite(i)-equi_pts_hermite(j)));
        end
    end
    
    f_k_i = f_prime_equi_pts_hemrite(i)*(grid_pts-equi_pts_hermite(i))*(l_i.^(2));
end

function [xx, ls, ws] = barycentric_interp(nodes, f_nodes)
    
end