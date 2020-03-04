clear
close all
clc
Q = 1/2*eye(2)*4;
R = 1/2*eye(1);

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

[K, P, eigenvalues] = dlqr(A, B, Q, R);

K;
eigenvalues;
poles = [0.9+0.1j; 0.9 - 0.1j];

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

t = [1:51];
figure(1);
plot(t, x(:,1), '-black');
hold on;
plot(t, x(:,2), '-black');
hold on;
plot(t, x_est(:,1), '-blue');
hold on;
plot(t, x_est(:,2), '-blue');


phi = [A-B*K B*K;
    zeros(2) A-Kf*c.']

eig(phi)
eig(A-B*K)
eig(A-Kf*c.')