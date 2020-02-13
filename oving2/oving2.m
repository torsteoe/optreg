c = -[100; 75; 55];
Aeq = [3, 2, 1; 2, 2, 3];
beq = [7200; 6001];
lb = zeros(3,1);
options = optimset('Algorithm', 'dual-simplex');
[x, fval, exitflag, output, lambda] = linprog(c, [],[], Aeq, beq, lb, [],[], options);
disp(x);
disp(fval);