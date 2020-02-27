function [x] = generateX()
%GENERATEX Summary of this function goes here
%   Detailed explanation goes here
    N=30;
    x0 = [0;0;1];
    x_real = zeros(N,1);
    
    
    x_real(1:3) = A*x0 + B*u(1);
    for i = 2:30
        x_real(i*3-2:i*3) = A*x_real(i*3-5:i*3-2) + B*u(i-1);
    end

end