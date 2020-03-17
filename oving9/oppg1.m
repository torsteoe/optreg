%% 1a) 1.2;1.2
close all;
clear;
clc;

x0 = [1.2; 1.2];


syms x1 x2 real;
x = [x1, x2];

f = @(x1, x2) 100*(x2-x1^2)^2 + (1-x1)^2;
xk = x0;
grad_f(x1, x2) = gradient(f, [x1,x2]);
hessian_f(x1,x2) = hessian(f, [x1, x2]);

p_k = newton_step(hessian_f, grad_f);
N = 20;
x_easy = zeros(2,N);
for k=1:N
    x_easy(:,k) = xk;
    step = double(subs(p_k, [x1, x2], xk'));
    alpha_1 = backtracking_line_search(f, simplify(grad_f), xk,step);
    xk = xk + alpha_1*step;
end

T = 1:N;


%% 1a) -1.2;1



x0 = [-1.2; 1];


syms x1 x2 real;
x = [x1, x2];

f = @(x1, x2) 100*(x2-x1^2)^2 + (1-x1)^2;
xk = x0;
grad_f(x1, x2) = gradient(f, [x1,x2]);
hessian_f(x1,x2) = hessian(f, [x1, x2]);

p_k = newton_step(hessian_f, grad_f);
N = 20;
x_difficult = zeros(2,N);
for k=1:N
    x_difficult(:,k) = xk;
    step = double(subs(p_k, [x1, x2], xk'));
    alpha_1 = backtracking_line_search(f, simplify(grad_f), xk,step);
    xk = xk + alpha_1*step;
end

T = 1:N;


%% 1b) BFGS

%Algorithm 6.1

x0 = [5; 5];
rho = 0.9;
dim = size(x0, 1);
k = 1;
H0 = eye(dim); %inverse hessian approximation
inverse_hessian = H0;
epsilon = 0.01;
alpha_k = 1;
xk = x0;
x_bfgs = x0;
grad_value_k = double(subs(grad_f, [x1, x2], xk'));
while norm(grad_value_k,2) > epsilon
    step = -inverse_hessian*grad_value_k;
    alpha_1 = backtracking_line_search(f, simplify(grad_f), xk,step);
    x_kp1 = xk + alpha_1*step;
    sk = x_kp1-xk;
    grad_value_k = double(subs(grad_f, [x1, x2], xk'));
    grad_value_kp1 = double(subs(grad_f, [x1, x2], x_kp1'));

    yk =  grad_value_kp1-grad_value_k;
    
    %compute next hessian
    inverse_hessian = inv_hessian(sk, yk, inverse_hessian);
    
    k = k+1
    xk = x_kp1
    grad_value_k = grad_value_kp1;
    x_bfgs = [x_bfgs, xk];
end

figure();
subplot(311);
plot(T, x_easy);
subplot(312);
plot(T, x_difficult);
subplot(313);
T_bfgs = 1:k;
plot(T_bfgs, x_bfgs);