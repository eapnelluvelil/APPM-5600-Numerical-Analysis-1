clc;
clear;

f = @(x) ones(size(x));
g = @(x) x;
h = @(x) x.^(5);


w = (1/6)*[1; 5; 5; 1];
xj = [-1; -1/sqrt(5); 1/sqrt(5); 1];

a = -1;
b = 1;

disp(integral(f, a, b));
disp(dot(f(xj), w));

disp(integral(g, a, b));
disp(dot(g(xj), w));

disp(integral(h, a, b));
disp(dot(h(xj), w));