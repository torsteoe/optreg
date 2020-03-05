
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
c = eye(2);
[K, P, eigenvalues] = dlqr(A, B, Q, R);


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
Beq(1:mx) = A*x0_hat;

G_with_p = eye(N*(mx+mu));

G_with_p(1:N*mx, 1:N*mx) = kron(eye(N), Q);
G_with_p(N*mx+1:N*mx+N*mu, N*mx+1:N*mx+N*mu) = kron(eye(N), R);
G_with_p((N-1)*(mx):N*mx-1, (N-1)*mx:N*mx-1) = G_with_p((N-1)*(mx):N*mx-1, (N-1)*mx:N*mx-1) + P;


LB = -Inf*ones(N*(mx+mu),1);
LB(N*mx+1:N*(mx+mu)) = -4;
UB = -LB;

options = optimset('Display', 'off', 'Diagnostics', 'Off', 'LargeScale', 'Off', 'Algorithm', 'interior-point-convex');
z = quadprog(G_with_p, [], [],[], Aeq, Beq, LB, UB);

% x_optimal_1 = z(1:2:N*mx);
% x_optimal_2 = z(2:2:N*mx);
% 
% 
%  plot(t, x_optimal_1, '-black');
%  hold on;
%  plot(t, x_optimal_2, '-black');
%  hold on;
poles = [0.5 + 0.03j; 0.5-0.03j];


Kf = place(A', c.', poles);


timesteps = 50;

current_state = x0;
xs_with_P = zeros(N, 2);

current_Beq = zeros(N*mx, 1);
current_Beq(1:mx) = A*x0;
inputs_with_p = zeros(timesteps*mu, 1);
for k = 1:timesteps
    
      xs_with_P(k,:) = current_state;

     current_Beq(1:mx) = A*current_state;
     [z,fval,exitflag,output,lambda] = quadprog(G_with_p, [], [],[], Aeq, current_Beq, LB, UB, current_state,[], options);
     current_input = z(N*mx+1);
     current_state = (A)*current_state + B*current_input;
     inputs_with_p(k) = current_input;

end



 t =1:timesteps;

subplot(211);
plot(t, xs_with_P(:,1), '-black');
hold on;
plot(t, xs_with_P(:,2), '-black');
hold on;


current_state = x0;
current_state_est = x0_hat;
current_input = 0;
xs = zeros(N, 2);
xs_est = zeros(N, 2);
G = eye(N*(mx+mu));
G(1:N*mx, 1:N*mx) = kron(eye(N), Q);
G(N*mx+1:N*mx+N*mu, N*mx+1:N*mx+N*mu) = kron(eye(N), R);
current_Beq = zeros(N*mx, 1);
current_Beq(1:mx) = A*x0;
timesteps = 50;
inputs = zeros(timesteps, 1);
for k = 1:timesteps
    xs(k,:) = current_state;

     current_Beq(1:mx) = A*current_state;
     [z,fval,exitflag,output,lambda] = quadprog(G, [], [],[], Aeq, current_Beq, LB, UB, current_state,[], options);
     current_input = z(N*mx+1);
     current_state = (A)*current_state + B*current_input;
     inputs(k) = current_input;

    
end

subplot(212);
plot(t, inputs_with_p, 'yellow');
hold on
plot(t, inputs, 'red');


subplot(211);
plot(t, xs(:,1), '-red');
hold on;
plot(t, xs(:,2), '-red');





