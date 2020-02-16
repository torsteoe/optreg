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

x_star = quadprog(G, [], [],[], Aeq, beq);

u = x_star(91:120);
x = x_star(1:90);
mx = 3;

%Removed x_0 
x1 = [x_star(1:mx:N*mx)];              % State x1 from solution
x2 = [x_star(2:mx:N*mx)];              % State x2 from solution
x3 = [x_star(3:mx:N*mx)]; 
delta_t = 0.25;
t = 0:delta_t:delta_t*(length(u)-1);

figure(2)
subplot(511)
stairs(t,u),grid
ylabel('u')
subplot(512)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('x1')
subplot(513)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('x2')
subplot(514)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('x3')

