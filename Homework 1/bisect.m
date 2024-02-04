%{
    Name: bisect
    Inputs: a - the left endpoint of the interval
            b - the right endpoint of the interval
            f - a continuous function defined on [a, b] such that f(a)f(b)
            < 0
            tol - the user-specified tolerance
                  When the approximate root, c, is such that |b-c| < tol,
                  bisect stops performing additional iterations
    Outputs: c - an approximate root of f on the interval [a, b]
             iters - the number of iterations it took for bisect to
             converge, which occurs either
                 - when f(c) == 0, to machine precision, or
                 - when |b-c| < tol
%}
function [c, iters] = bisect(a, b, f, tol)
    iters = 1;
    % Compute first approximation to root
    c = (a+b)/2;
    % While the interval of interest is large enough and
    % f evaluated at the current iterate is not 0,
    % perform bisection
    while ((f(c) ~= 0) && (abs(b-c) >= tol))
        iters = iters + 1;
        if (f(a)*f(c) < 0)
            b = c;
        else
            a = c;
        end
        c = (a+b)/2; 
    end
end