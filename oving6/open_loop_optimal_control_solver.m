clear;
clc;
close all;
Q = [0 0 0;
    0 0 0;
    0 0 2];
r = 1;
R = 2*r;
N = 30;
A = [0 0 0;
    0 0 1;
    0.1 -0.79 1.78];
b = [1;0;0.1];
x_0 = [0;0;1];

q = diag_repeat(Q, N);
p = diag_repeat(R, N);
G = blkdiag(q,p);
c = zeros(120,1);

a1 = zeros(N*3, N*3);
b1	= diag_repeat(-b,N);
a1(1:3, 1:3) = eye(3);
for k = 2:N
    a1((k*3-2):(k*3-2+2), (k*3-2):(k*3-2+2)) = eye(3);
    a1((3*k-2):((3*k-2)+2), (3*k-5):(3*k-5)+2) = -A;
end

Aeq = [a1 b1];
beq = zeros(N*3, 1);
beq(1:3,1) = A*x_0;


mx = 3;
LB = [zeros(90, 1); -ones(30, 1)];
UB = [zeros(90, 1); ones(30, 1)];
for k = 1:90
    LB(k) = -Inf;
    UB(k) = Inf;
end
[x_star2, Fval, exitflag, outp] = quadprog(G, [],[],[], Aeq, beq, LB, UB);
u_quad = x_star2(91:120,1);

%Removed x_0 
delta_t = 0.25;
t = 0:delta_t:delta_t*(length(u_quad)-1);

x1_quad = [x_star2(1:mx:N*mx)];              % State x1 from solution
x2_quad = [x_star2(2:mx:N*mx)];              % State x2 from solution
x3_quad = [x_star2(3:mx:N*mx)]; 





figure(2)
subplot(211)
plot(t, x3_quad, '+-'),grid
ylabel('x3')
subplot(212)
plot(t, u_quad, '+'),grid




