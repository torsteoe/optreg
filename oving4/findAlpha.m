A = [1 -2;
    -1 -2;
    -1 2;
    1 0;
    0 1];

b = [-2;-6;-2; 0; 0];

p = [-1.2;2.4];
x = [2.2;0.1];

above = (b-A*x);
below = A*p;
disp(above);
disp(below);

x3 = x + 2/3*p;
disp(x3);

