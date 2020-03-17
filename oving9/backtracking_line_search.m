function [alpha_k] = backtracking_line_search(f, grad_f, x_k, p_k, alpha_bar, rho, c)
%BACKTRACKING_LINE_SEARCH Summary of this function goes here
%   Detailed explanation goes here


if nargin < 5
    alpha_bar = 1;
end
if nargin < 6
    rho = 0.9;
end
if nargin < 7
    c=1e-4;
end
alpha_k = alpha_bar;
cell = num2cell(x_k);
temp_step = x_k+alpha_k*p_k ;
double(c*alpha_k*grad_f(cell{:}));
while f(temp_step(1), temp_step(2))  > f(x_k(1), x_k(2)) +double(c*alpha_k*grad_f(cell{:})).'*p_k
    
    alpha_k = rho*alpha_k;
    temp_step = x_k+alpha_k*p_k ;
end

end

