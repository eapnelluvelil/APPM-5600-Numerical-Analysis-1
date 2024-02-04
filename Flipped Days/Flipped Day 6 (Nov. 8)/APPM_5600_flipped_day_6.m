% Problem 1
a = -1;
b = 1;
num_pts = 1000;
mesh = linspace(a, b, num_pts)';

n_final = 3;

Legendre_polys = zeros(num_pts, n_final);
Legendre_polys(:, 1) = ones(num_pts, 1);
Legendre_polys(:, 2) = mesh;

for n = 3:n_final
    Legendre_polys(:, n) = ((2*n+1)/(n+1))*mesh.*Legendre_polys(:, n-1) - (n/(n+1))*Legendre_polys(:, n-2);

    figure(1);
    plot(mesh, Legendre_polys(:, n), 'LineWidth', 2, 'DisplayName', string(n-1));
    hold on;
end

figure(1);
legend;
ylim([-1, 1]);

