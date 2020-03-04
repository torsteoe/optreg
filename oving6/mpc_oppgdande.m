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
G_without_ib = blkdiag(q,p);
c = zeros(120,1);

a1 = zeros(N*3, N*3);
b1	= diag_repeat(-b,N);
a1(1:3, 1:3) = eye(3);
for k = 2:N
    a1((k*3-2):(k*3-2+2), (k*3-2):(k*3-2+2)) = eye(3);
    a1((3*k-2):((3*k-2)+2), (3*k-5):(3*k-5)+2) = -A;
end

Aeq_without_ib = [a1 b1];
beq_without_ib = zeros(N*3, 1);
beq_without_ib(1:3,1) = A*x_0;
beq_with_ib = beq_without_ib;

mu = 3;
b1 = zeros(N*3, 6);
b1(1:mu) = b;
b1(mu+1:mu+mu,2) = b;
b1(mu*2+1:mu*2+2*mu,3) = [b;b];
b1(mu*4+1:mu*4+4*mu,4) = kron(ones(4,1), b);
b1(mu*8+1:mu*8+8*mu,5) = kron(ones(8,1), b);
b1(mu*16+1:mu*16+14*mu,6) = kron(ones(14,1), b);

Aeq_with_ib = [a1 -b1];
p_with_ib = diag_repeat(R*5, 6);
G_with_ib = blkdiag(q,p_with_ib);
mx = 3;
LB = [zeros(90, 1); -ones(30, 1)];
UB = [zeros(90, 1); ones(30, 1)];
LB_with_ib = [zeros(90, 1); -ones(6, 1)];
UB_with_ib = [zeros(90, 1); ones(6, 1)];
for k = 1:90
    LB(k) = -Inf;
    LB_with_ib(k) = -Inf;
    UB(k) = Inf;
    UB_with_ib(k) = Inf;
end

u_quad_mpc = zeros(N,1);
u_quad_mpc_with_ib = zeros(N,1);

state_without_ib = x_0;
states_without_ib = zeros(N*3,1);
states_without_ib(1:3) = x_0;
state_with_ib = state_without_ib;
states_with_ib = states_without_ib;
for iteration = 1:N
    
%mpc without input blocking    
% [x_star2, Fval, exitflag, outp] = quadprog(G_without_ib, [],[],[], Aeq_without_ib, beq_without_ib, LB, UB);
% first_u = x_star2(91);
% u_quad_mpc(iteration) = first_u;
% state_without_ib = A*state_without_ib + b*first_u;
% beq_without_ib(1:3,1) = A*state_without_ib;
% states_without_ib(iteration*3+1:iteration*3+3) = state_without_ib;

%mpc with input blocking
[x_star2_ib, Fval, exitflag, outp] = quadprog(G_with_ib, [],[],[], Aeq_with_ib, beq_with_ib, LB_with_ib, UB_with_ib);
first_u = x_star2_ib(91);
u_quad_mpc_with_ib(iteration) = first_u;
state_with_ib = A*state_with_ib + b*first_u;
beq_with_ib = zeros(N*3, 1);
beq_with_ib(1:3,1) = A*state_with_ib;
states_with_ib(iteration*3+1:iteration*3+3) = state_with_ib;



end
%input blocking scheme without mpc
beq_with_ib = zeros(N*3, 1);
beq_with_ib(1:3,1) = A*x_0;
[x_star2, Fval, exitflag, outp] = quadprog(G_with_ib, [],[],[], Aeq_with_ib, beq_with_ib, LB_with_ib, UB_with_ib);
%u_quad = x_star2(91:N*3+blocks);
u_quad = zeros(N, 1);
for k = 1:N
    i = floor(k/6);
    u_elem = x_star2(N*3+1+i);
    u_quad(k) = u_elem;
end

% State x2 from solution
x3_quad = [x_0(3);x_star2(3:mx:N*mx-3)]; 








%Removed x_0 
delta_t = 0.25;
t = 0:delta_t:delta_t*(length(u_quad_mpc_with_ib)-1);





figure(2)
subplot(211)
plot(t, states_with_ib(3:3:N*3), 'o', t, x3_quad, '-'),grid
ylabel('x3')
subplot(212)
plot(t, u_quad_mpc_with_ib, 'o', t, u_quad, '-'),grid




