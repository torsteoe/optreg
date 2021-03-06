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
mat = [G -Aeq';
        Aeq zeros(N*3)];
c_b = [-c;beq];

A_ineq = [zeros(30, 30*3) eye(30);
          zeros(30, 30*3) -eye(30)];
b_ineq = [ones(30,1); ones(30,1)];
answer = mat\c_b;

x_star = answer(1:N*3+N,1);
lambda_star = answer(N*3+30:end,1);

u = x_star(91:120);
x = x_star(1:90);
mx = 3;
[x_star2, Fval, exitflag, outp] = quadprog(G, [],[],[], Aeq, beq);
u_quad = x_star2(91:120);

%Removed x_0 
x1 = [x_star(1:mx:N*mx)];              % State x1 from solution
x2 = [x_star(2:mx:N*mx)];              % State x2 from solution
x3 = [x_star(3:mx:N*mx)]; 
delta_t = 0.25;
t = 0:delta_t:delta_t*(length(u)-1);

x1_quad = [x_star2(1:mx:N*mx)];              % State x1 from solution
x2_quad = [x_star2(2:mx:N*mx)];              % State x2 from solution
x3_quad = [x_star2(3:mx:N*mx)]; 
u_combined = [u u_quad];

x_mpc = zeros(N,1);
x_open_loop = zeros(N, 1);
    
x_mpc(1:3) = A*x_0 + b*u(1);
u_mpc = zeros(30,1);
u_mpc(1) = u_quad(1);


A_plant = [0 0 0;
           0 0 1;
           0.1 -0.855 1.85];
b_plant = [1;0;0];

for i = 2:30
    
    x_this_timestep_no_mpc = A_plant*x_open_loop(i*3-5:i*3-3) + b_plant*u_quad(i-1);
    x_open_loop(i*3-2:i*3) = x_this_timestep_no_mpc;

    
    x_this_timestep_mpc_controlled =  A_plant*x_mpc(i*3-5:i*3-3) + b_plant*u_mpc(i-1);    
    beq(1:3,1) = A*x_this_timestep_mpc_controlled;
    z_mpc = quadprog(G,[],[],[], Aeq, beq);
    u_mpc(i) = z_mpc(N*mx+1);
    x_mpc(i*3-2:i*3) = x_this_timestep_mpc_controlled;
end

x_real_3_mpc = [x_mpc(3:mx:N*mx)];
x_real_3_open_loop = [x_open_loop(3:mx:N*mx)];

figure(2)
subplot(511)
stairs(t,u_combined),grid
ylabel('u')
subplot(512)
plot(t,x1,'m', t, x1_quad, '+'),grid
ylabel('x1')
subplot(513)
plot(t,x2,'m', t, x2_quad, '+'),grid
ylabel('x2')
subplot(514)
plot(t,x3,'m', t, x3_quad, '+'),grid
ylabel('x3')
subplot(515)
plot(t,x_real_3_mpc,'m', t, x_real_3_open_loop, 'o'),grid
ylabel('x_{real}')


