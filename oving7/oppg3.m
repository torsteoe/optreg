
clear
close all
clc
Q = eye(2)*4;
R = eye(1);

k1 = 1;
k2 = 1;
k3 = 1;
T = 0.1;

A = [1 T;
    -k2*T 1-k1*T];

B = [0;k3*T];
x0 = [5;1];
x0_hat = [6;0];
c = [1 0];

[K, P, eigenvalues] = dlqr(A, B, Q/2, R/2);

K;
eigenvalues;
poles = [0.5+0.03j; 0.5 - 0.03j];

mx = 2;
mu = 1;
N = 10;

%Make Aeq

A1 = eye(N*mx);

for k=1:N-1
    A1(k*mx+1:k*mx+mx,k*mx-1:k*mx) = -A;
end

B1 = kron(eye(N), -B);
Aeq = [A1 B1];

Beq = zeros(N*mx, 1);
Beq(1:mx) = A*x0;

G = eye(N*(mx+mu));

G(1:N*mx, 1:N*mx) = kron(eye(N), Q);
G(N*mx+1:N*mx+N*mu, N*mx+1:N*mx+N*mu) = kron(eye(N), R);

z0 = zeros(N*(mx+mu), 1);
z0(1:mx) = x0;
LB = -Inf*ones(N*(mx+mu),1);
LB(N*mx+1:N*(mx+mu)) = -4;
UB = -LB;

options = optimset('Display', 'off', 'Diagnostics', 'Off', 'LargeScale', 'Off', 'Algorithm', 'interior-point-convex');
z = quadprog(G, [], [],[], Aeq, Beq, LB, UB, x0_hat);


x_optimal_1 = z(1:2:N*mx);
x_optimal_2 = z(2:2:N*mx);


% plot(t, x_optimal_1, '-black');
% hold on;
% plot(t, x_optimal_2, '-black');
% hold on;

Kf = place(A, c.', poles);

current_state = x0;
current_state_est = x0_hat;
current_input = 0;
xs = zeros(N, 2);
xs_est = zeros(N, 2);

current_Beq = zeros(N*mx, 1);
current_Beq(1:mx) = A*x0;
for k = 1:50
    xs(k,:) = current_state;
    xs_est(k,:) = current_state_est;
    if mod(k, 5) == 0 
        current_Beq(1:mx) = A*current_state_est;
        z = quadprog(G, [], [],[], Aeq, current_Beq, LB, UB, current_state_est, options);
        current_input = z(N*mx+1);
    end
     y = c*current_state;
     
     y_hat = c*current_state_est;
     current_state_est = (A)*(current_state_est) + B*current_input +Kf.'*(y-y_hat);
     current_state = (A)*current_state + B*current_input;

     

%     
end

timesteps = 50;
t =1:timesteps;
plot(t, xs(:,1), '-black');
hold on;
plot(t, xs(:,2), '-black');
hold on;
plot(t, xs_est(:,1), '-blue');
hold on;
plot(t, xs_est(:,2), '-blue');




Kf = place(A, c.', poles);
x = zeros(51, 2);
x(1,:) = x0;
x_est = zeros(51, 2);
x_est(1,:) = x0_hat;


for k = 1:50
    x(k+1,:) = (A-B*K)*x(k,:).';
    y = c*(x(k,:)).';
    y_hat = c*(x_est(k,:)).';
    x_est(k+1,:) = (A-B*K)*(x_est(k,:)).'+Kf.'*(y-y_hat);
end
% 
% t = [1:51];
% figure(1);
% plot(t, x(:,1), '-black');
% hold on;
% plot(t, x(:,2), '-black');
% hold on;
% plot(t, x_est(:,1), '-blue');
% hold on;
% plot(t, x_est(:,2), '-blue');





phi = [A-B*K B*K;
    zeros(2) A-Kf*c.'];

eig(phi);
eig(A-B*K);
eig(A-Kf*c.');
