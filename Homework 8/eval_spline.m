%{
    Function: eval_spline
    Inputs:
        a - left endpoint of interval (x(1) = a)
        b - right endpoint of interval (x(end) = b)
        c - vector of coefficients for B-splines
        x - coarse grid containing sub-intervals [x_j, x_{j+1}]
        t_plot_m - number of fine data points to use when plotting
    Outputs:
        t_plot - vector of fine data points (length(t_plot) = t_plot_m)
        f_t_plot - vector of spline evaluated at t_plot
%}
function [t_plot, f_t_plot] = eval_spline(a, b, c, x, t_plot_m)
    x_h = abs(x(2) - x(1));
    t_plot = linspace(a, b, t_plot_m)';
    f_t_plot = zeros(size(t_plot));

    for i = 1:t_plot_m
        for j = 1:length(x)
            f_t_plot(i) = f_t_plot(i) + c(j)*eval_B_spline(x(j), t_plot(i), x_h);
        end
    end
end