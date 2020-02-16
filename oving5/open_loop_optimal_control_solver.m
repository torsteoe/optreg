Q = [0 0 0;
    0 0 0;
    0 0 2];
r = 1;
R = 2*r;
N = 30;
q = diag_repeat(Q, N);
p = diag_repeat(R, N);
G = blkdiag(q,p);