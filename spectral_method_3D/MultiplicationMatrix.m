function M = MultiplicationMatrix(v,lambda)
% returns the multiplication matrix in C^{lambda} basis for the function
% described by the chebyshev coefficients in v
n = length(v);
M = zeros(n,n);
if lambda == 0
    v1 = v; v1(1) = 2*v1(1);
    M = toeplitz(v1);
    N = hankel(v);
    M(2:end,:) = M(2:end,:) + N(2:end,:);
    M = M./2;    
   return 
end
if lambda == 1
    v1 = v; v1(1) = 2*v1(1);
    M = toeplitz(v1);
    M(1:end-2, 1:end-2) = M(1:end-2, 1:end-2) - hankel(v(3:end));
    M = M./2;
   return 
end
d1 = (1:(n-1))./(2*(lambda:lambda+n-2));
d2 = (2*lambda+(0:n-2))./(2*(lambda+1:lambda+n-1));
N = diag(d1,-1)+diag(d2,1); 
for k =  0:lambda-1 
v = getSUltra(n,k)*v; 
end
A = cell(n);
A{1} = eye(n);
A{2} = 2*lambda*N;
for i = 1:n-2
    A{i+2} =  2*(i+lambda)/(i+1)*N*A{i+1}-(i+2*lambda -1)/(i+1)*A{i};
end
for i = 1:n
    M = M + v(i).*A{i};
end
end