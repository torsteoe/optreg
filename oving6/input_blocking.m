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
blocks = 6;

q_open_loop = diag_repeat(Q, N);
p_open_loop = diag_repeat(R*5, blocks);
G = blkdiag(q_open_loop,p);
c = zeros(120,1);

a1 = zeros(N*3, N*3);
b1	= zeros(N*3, blocks);
num_states_same_input = 3*5;
for k = 1:blocks
    b1((k-1)*(num_states_same_input)+1:(k-1)*(num_states_same_input)+num_states_same_input, k) = -kron(ones(5, 1), b);
    
end

a1(1:3, 1:3) = eye(3);
for k = 2:N
    a1((k*3-2):(k*3-2+2), (k*3-2):(k*3-2+2)) = eye(3);
    a1((3*k-2):((3*k-2)+2), (3*k-5):(3*k-5)+2) = -A;
end

Aeq = [a1 b1];
beq = zeros(N*3, 1);
beq(1:3,1) = A*x_0;


mx = 3;
LB = [ones(90, 1)*-Inf; -ones(6, 1)];
UB = -LB;
for k = 1:90
    LB(k) = -Inf;
    UB(k) = Inf;
end
[x_star2, Fval, exitflag, outp] = quadprog(G, [],[],[], Aeq, beq, LB, UB);
%u_quad = x_star2(91:N*3+blocks);
u_quad = zeros(N, 1);
for k = 1:N
    i = floor(k/6);
    u_elem = x_star2(N*3+1+i);
    u_quad(k) = u_elem;
end
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




