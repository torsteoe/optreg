clf;
clc;
clear;


%Creating contour, objective function:
linx = linspace(-5,15);
liny = linspace(-5,15);
[x1, x2] = meshgrid(linx,liny);


c = [-3 -2 0 0]';
A = [2 1 1 0;
    1 3 0 1];
b = [8; 15];
obj = -3*x1 -2*x2;
figure(1);
hold on
xlabel("x_1"); ylabel("x_2");
contour(x1, x2, obj);

inequality_1 = 2*x1 + x2 <= 8;
inequality_2 = x1 + 3*x2 <= 15;
inequality_3 = x1>=0;
inequality_4 = x2>=0;

[c1, h1] = contour(x1, x2, inequality_1);
h1.LineColor = "black";
[c2, h2] =contour(x1, x2, inequality_2);
h2.LineColor = "black";

[c3, h3] =contour(x1, x2, inequality_3);
h3.LineColor = "black";

[c4, h4] =contour(x1, x2, inequality_4);
h4.LineColor = "black";


%2c)

x0 = [0 0 8 15]'; % The starting point must also have four variables!

% Solve problem and print a report: (Call simplex(c,A,b,x0,[]) if you do 
% not want the report printed).
[x, fval, iterates] = simplex(c,A,b,x0,'report');

% Extract iterates in the space of x_1 and x_2:
iter_x1_x2 = iterates(1:2, :);

% Extract iterates as individual vectors (containing x_1 and x_2):
iter_1 = iter_x1_x2(:,1);
iter_2 = iter_x1_x2(:,2);
iter_3 = iter_x1_x2(:,3);

plot(iter_1(1), iter_1(2), 'r-o', 'Markersize', 12);
text(iter_1(1)+0.3, iter_1(2)+0.3, "step 1")
plot(iter_2(1), iter_2(2), 'r-o', 'Markersize', 12);
text(iter_2(1)+0.3, iter_2(2)+0.3, "step 2")
plot(iter_3(1), iter_3(2), 'r-o', 'Markersize', 12);
text(iter_3(1)+0.3, iter_3(2)+0.3, "step 3")

plot(1.8, 4.4, 'yellow-x');
