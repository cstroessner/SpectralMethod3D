function A = getForwardOperatorMatrix(opLCP)
% returns matrix A such that LuCoeffs = reshape(A*uCoeffs(:),n);
r = size(opLCP,2);
n = [size(opLCP{1,1},2),size(opLCP{2,1},2),size(opLCP{3,1},2)];
A = sparse(n(1)*n(2)*n(3),n(1)*n(2)*n(3));
for i = 1:r
    A = A + kron(opLCP{3,i},kron(opLCP{2,i},opLCP{1,i}));
end
end