function [C, gamma] = MaxPH(A, B, alpha, betta)
n = length(A);
m = length(B);
a = -A*ones(n, 1);
b = -B*ones(m, 1);
C = [KronSum(A,B),   kron(eye(n),b),  kron(a, eye(m)); 
    zeros(n, n*m),     A,               zeros(n, m);
    zeros(m,n*m),     zeros(m, n),   B] ;
alpha0 = 1-sum(alpha);
betta0 = 1 - sum(betta);
gamma = [kron(alpha, betta), betta0*alpha, alpha0*betta];