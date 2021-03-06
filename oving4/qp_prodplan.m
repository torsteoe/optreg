% ********************************************************
% *                                                      *
% *      Optimalisering og regulering 					 *
% * 		?ving 3 Oppgave   V?r 2003					 *
% *                                                      *
% *      Bjarne Foss 1996                                *
% *                                                      *
% * qp_prodplan.m                                        *
% *                                                      *
% * m-file for calculating QP solution.                  *
% *                                                      *
% * Oppdated 10/1-2001 by Geir Stian Landsverk           *
% *                                                      *
% * Verified to work with MATLAB R2015a,                 *
% *    Andreas L. Fl?ten                                 *
% *                                                      *
% * Verified to work with MATLAB R2018b,                 *
% *    Joakim R. Andersen                                *
% *                                                      *
% ********************************************************
clf;
clc;
clear;

global XIT; % Storing iterations in a global variable
global IT;  % Storing number of iterations in a global variable
global x0;  % Used in qp1.m (line 216).

IT=1; XIT=[];

x0 = [0 0]'; % Initial value P1
vlb = [0 0]; % Lower bound on x
vub = [];    % Upper bound on x

% min 0.5*x'*G*x + x'*c
%  x 
%
% s.t. A*x <= b 

% Quadratic objective (MODIFY THESE)
G = [0.8 0;
     0 0.4]; % Remember the factor 1/2 in the objective
c = [-3 -2];

% Linear constraints (MODIFY THESE)
A = [2 1;
    1 3];
b = [8; 15];

options = optimset('LargeScale','Off');
[x,lambda] = quadprog1(G,c,A,b,[],[],vlb,vub,x0,options);

disp('Iteration sequence:')
disp(XIT');
disp('Solution:')
disp(x);
A*x-b
disp(lambda);
hold on
linx = linspace(-5,10);
liny = linspace(-5,10);
[x1,x2] = meshgrid(linx, liny);

obj = -(x1*3-x1.^2*0.4)-(x2*2-x2.^2*0.2);
contour(x1, x2, obj, 20);

ineq1 = x1>=0;
ineq2 = x2>=0;
ineq3 = 2*x1 + x2 <=8;
ineq4 = 1*x1 + 3*x2 <=15;


[c1, h1] = contour(x1, x2, ineq1);
h1.LineColor = "black";
[c2, h2] =contour(x1, x2, ineq2);
h2.LineColor = "black";

[c3, h3] =contour(x1, x2, ineq3);
h3.LineColor = "black";

[c4, h4] =contour(x1, x2, ineq4);
h4.LineColor = "black";
for i=1:3
    plot(XIT(i, 1), XIT(i,2), 'yellow-x')
end





