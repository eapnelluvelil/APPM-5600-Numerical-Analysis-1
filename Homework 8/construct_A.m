%{
%}
function [c, x] = construct_A(a, b, num_x_intervals, m, f, noise)
    x   = linspace(a, b, num_x_intervals + 1)';
    x_h = abs(x(2)-x(1));

    % Include one equally spaced point before x_0 and one equally spaced
    % point after x_n
    x = [x(1) - x_h; ...
         x; ...
         x(end) + x_h];

    t = linspace(a, b , m)';

    A = zeros(m, num_x_intervals + 3);

    for i = 1:m
        for j = 1:length(x)
            A(i, j) = eval_B_spline(x(j), t(i), x_h);
        end
    end

    figure;
    spy(A'*A);
    title("Spy plot of A^{T}A")

    rhs = f(t) + noise;
    c = A\rhs;
end