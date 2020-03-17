function [p_k] = newton_step(hessian, gradient)
%NEWTON_STEP Summary of this function goes here
%   Detailed explanation goes here
p_k = -hessian\gradient;

end

