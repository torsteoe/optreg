function [next_hessian] = inv_hessian(s, y, prev_hessian)
%INV_HESSIAN Compute inverse hessian from previous inverse hessian
%6.17 in Numerical optimization
rho = 1/(y.'*s);
I = eye(size(prev_hessian, 2));
next_hessian = (I-rho*s*y.')*prev_hessian*(I-rho*y*s.') + rho*(s*s.');

end

