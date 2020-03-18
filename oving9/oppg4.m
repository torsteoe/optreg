%% Testing Nelder Mead
close all;
clear;
clc;
x0 = [-2;2];
f1 = figure();
[x, fval, iter] = nelder_mead(x0, ' ');
legend('[-2,2]');
x0 = [-1.2; 1];
f2 = figure();
[x, fval, iter] = nelder_mead(x0, ' ');
legend('[-1.2, 1]');

x0 = [5;5];
f3 = figure();
[x, fval, iter] = nelder_mead(x0, ' ');
legend('[5,5]');

x0 = [1;1];
f4 = figure();
[x, fval, iter] = nelder_mead(x0, ' ');
legend('[1,1]');

movegui(f1,'west');
movegui(f2,'north');
movegui(f3,'east');
movegui(f4,'south');

%% 1c
x0 = [1.2;1.2];
f1 = figure();
[x, fval, iter] = nelder_mead(x0, ' ');
legend('[1.2,1.2]');
f2 = figure();
plot_iter_rosenbrock(iter);
legend('[1.2,1.2]');

x0 = [-1.2; 1];
f3 = figure();
[x, fval, iter] = nelder_mead(x0, ' ');
legend('[-1.2, 1]');


f4 = figure();
plot_iter_rosenbrock(iter(:, 1:size(iter, 2)-1));
legend('[-1.2, 1]');


movegui(f1,'west');
movegui(f2,'north');
movegui(f3,'east');
movegui(f4,'south');