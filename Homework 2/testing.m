a = 3;
b = 4;
c = 5;

f_coeffs = [1, -1*(a+b+c), (c*(a+b) + a*b), -a*b*c];
% f = @(x) polyval(f_coeffs, x);
f = @(x) (x-a).*(x-b).*(x-c);
f_prime_coeffs = [3, -2*(a+b+c), (c*(a+b) + a*b)];
f_prime = @(x) polyval(f_prime_coeffs, x);
% f_prime = @(x) 3*x.^(2) - 2*(a + b + c)*x + (c*(a+b) + a*b);

tol = 1e-13;
n_max = 100;
x_0 = -100000;

[r, iters, ~] = newton(x_0, f, f_prime, tol, n_max);
fprintf("Root computed via Newton's method: %0.16f (%d iterations, x_0 = %0.16f)\n", ...
        r, iters, x_0);
    
f_prime_roots = roots(f_prime_coeffs);
x_0 = f_prime_roots(2);
% [r, iters, ~] = newton(x_0, f, f_prime, tol, n_max);
% fprintf("Root computed via Newton's method: %0.16f (%d iterations, x_0 = %0.16f)\n", ...
%         r, iters, x_0);