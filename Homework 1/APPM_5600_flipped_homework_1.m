fprintf("Flipped homework 1\n");

main();

function main()
    funcs = {@(x) (x-1).*(x-3).*(x-5);
             @(x) ((x-1).^(2)).*(x-3);
             @(x) sin(x)
             };
    a_vals = [0; 0; 0];
    b_vals = [2.4; 2; 0.1];
    eps    = 1e-5;
    for i=1:length(funcs)
        c = bisect(a_vals(i), b_vals(i), funcs{i}, eps);
        fprintf("root r = %0.16f, f(r) = %0.16f\n", c, funcs{i}(c));
    end
end