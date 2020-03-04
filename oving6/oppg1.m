A = [1 0.5;
    0 1];
b = [0.125; 0.5] ;

Q = eye(2); % divided by two because dlqr
R = 1;

[K, P, eigenvalues] = dlqr(A, b, Q, R)