%% Flipped day 4
clear;
close all;

%% Problem 2
space = [0.01, 0.1];
x = @(x, y, h) 1 + h^(2)*(exp(y.*sqrt(x)) + 3*x.^(2));
y = @(x, y, h) 0.5 + h^(2)*tan(exp(x) + y.^(2));

x0 = [1; 0.5];
iters = 1;
for i=1:length(space)
    x1 = [x(x0(1), x0(2), space(i)); y(x0(1), x0(2), space(i))];
    while norm(x1-x0, Inf) >= 0.1*(space(i))^(4)
        x0 = x1;
        iters = iters + 1;
    end
    fprintf("Number of iterations required with h = %0.16f: %d\n", ...
        space(i), iters);
end

J = @(x, y, h) [h^(2)*((1/2)*(y./sqrt(x)).*exp(y.*sqrt(x)) + 6*x), h^(2)*sqrt(x).*exp(y.*sqrt(x)); ...
    h^(2)*(sec(exp(x) + y.^(2)).^(2).*exp(x)), h^(2)*sec(exp(x) + y.^(2)).^(2).*(2*y)];

for i=1:length(space)
    fprintf("");
end