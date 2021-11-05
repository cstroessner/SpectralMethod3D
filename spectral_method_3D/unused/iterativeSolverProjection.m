function [x,beta] = iterativeSolverProjection(Op,b)
% this uses GMRES as iterative solver (matrix multiplication is performed
% using the fast apply forward operator function)
k = 30;
m = 30000;
tol = 1e-15;

n = size(b);
if ~exist('x0','var')
    x0 = zeros(n);
end

i = 0;
beta = [];

while i < m
    if i+k > m
        k = m-i;
    end
    [x, beta_new] =opGMRES(Op,b,x0,k,tol);
    beta = [beta, beta_new];
    
    if beta(end) < tol * beta(1)
        i = i + length(beta)-1;
        fprintf('GMRES converged after %i steps.\n', i);
        break
    end
    x0 = x;
    i = i + k;    
end
if i >= m
    warning('GMRES did not converge %.2e.\n', beta(end))
end
end

function [x,beta] = opGMRES(Op,b,x0,k,tol)
% perform at most k steps of GMRES or until tol is reached
n = size(x0);
r = b-applyForwardOperator(Op,x0);
beta0 = norm(reshape(r,[prod(n),1]));
beta = beta0;
e1 = zeros(k+1,1);
e1(1) = beta0;
U = zeros(prod(n), k+1);
U(:,1) = r(:) / beta0;
H = zeros(k+1, k);
for j=1:k
    w = reshape(projectHomogDirichlet(applyForwardOperator(Op,reshape(U(:,j),n))),[prod(n),1]);
    H(1:j,j) = U(:,1:j)'*w;
    u_new = w - U(:,1:j)*H(1:j,j);
    H(j+1,j) = norm(u_new);
    u_new = u_new / H(j+1,j);
    U(:,j+1) = u_new;
    
    y = H(1:j+1,1:j) \ e1(1:j+1,1);
    beta = [beta, norm( e1(1:j+1,1) - H(1:j+1,1:j)*y )];
    
    if beta(end) < tol * beta0
        %fprintf('GMRES converged after %i steps.\n', j);
        break
    end
    if j == k
        %warning('GMRES did not converge %.2e.\n', beta(end))
    end
end
x = x0 + reshape(U(:,1:j)*y,n);
end