% Homework 9

clc;
clear;
close all;

n = 3;

S = diag(ones(2*n, 1), 1);
S(end, 1) = 1;

E = zeros(2*n+1);
for j = 1:(2*n+1)
    E(:, j) = exp((0:2*n)'*(2*pi*1i*j)/(2*n+1));
end

a = (1:(2*n+1))';

C = zeros(2*n+1);
for i=0:2*n
    C = C + a(i+1)*S^(i);
end

x = (1:(2*n+1))';

