function [val] = rosenbrock(x)
    x1 = x(1);
    x2 = x(2);
    val = 100*(x2-x1^2)^2 + (1-x1)^2;
end
