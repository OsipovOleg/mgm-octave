function C= KronSum(A, B)
n = length(A);
m = length(B);
C = kron(A, eye(m))+kron(eye(n), B);
