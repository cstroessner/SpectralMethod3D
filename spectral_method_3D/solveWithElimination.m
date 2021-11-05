function uCoeffs = solveWithElimination(opLCP,D,T1,F1,T2,F2,T3,F3,fastFlag)
% This function sets up all terms for the reduced system after substitution and calls
% solveSylvester to compute uCoeffs(k+1:end,k+1:end,k+1:end). The remaining terms are added based on
% the boundary conditions. In the end we return the solution of the PDE.

n = size(D); uCoeffs = zeros(n); k = 2; r = size(opLCP,2);

% compute the modified rhs for the reduced system
for j = 1:r
    D = D - tprod(F1, opLCP{1,j}(:,1:k), opLCP{2,j}, opLCP{3,j});
    D = D - tprod(F2, opLCP{1,j} - opLCP{1,j}(:,1:k)*T1, opLCP{2,j}(:,1:k), opLCP{3,j});
    D = D - tprod(F3, opLCP{1,j} - opLCP{1,j}(:,1:k)*T1, opLCP{2,j} - opLCP{2,j}(:,1:k)*T2 ,opLCP{3,j}(:,1:k));
end
D = D(1:end-k,1:end-k,1:end-k);

% modify the lhs and get the reduced system
for i = 1:r
opLCP{1,i} = opLCP{1,i} - sparse(opLCP{1,i}(:,1:k)*T1);
opLCP{1,i} = opLCP{1,i}(1:end-k,k+1:end);
opLCP{2,i} = opLCP{2,i} - sparse(opLCP{2,i}(:,1:k)*T2);
opLCP{2,i} = opLCP{2,i}(1:end-k,k+1:end);
opLCP{3,i} = opLCP{3,i} - sparse(opLCP{3,i}(:,1:k)*T3);
opLCP{3,i} = opLCP{3,i}(1:end-k,k+1:end);
end

% solve Sylvester equation
uCoeffs(k+1:end,k+1:end,k+1:end) = solveLinearEquation(opLCP,D,fastFlag);

% find the other blocks by using the boundary conditions
uCoeffs(1:2,k+1:end,k+1:end) = F1(:,k+1:end,k+1:end) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,k+1:end,k+1:end),[n(1)-k,(n(2)-k)*(n(3)-k)]),[2,n(2)-k,n(3)-k]);
uCoeffs(k+1:end,1:2,k+1:end) = F2(k+1:end,:,k+1:end) - permute(reshape(T2(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,k+1:end),[2,1,3]),[n(2)-k,(n(1)-k)*(n(3)-k)]),[2,n(1)-k,n(3)-k]),[2,1,3]);
uCoeffs(k+1:end,k+1:end,1:2) = F3(k+1:end,k+1:end,:) - permute(reshape(T3(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,k+1:end),[3,2,1]),[n(3)-k,(n(2)-k)*(n(1)-k)]),[2,n(2)-k,n(1)-k]),[3,2,1]);
uCoeffs(1:k,k+1:end,1:k) = F1(:,k+1:end,1:k) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,k+1:end,1:2),[n(1)-k,(n(2)-k)*k]),[2,n(2)-k,k]);
uCoeffs(1:k,1:k,k+1:end) = F1(:,1:k,k+1:end) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,1:2,k+1:end),[n(1)-k,(n(3)-k)*k]),[2,k,n(3)-k]);
uCoeffs(k+1:end,1:k,1:k) = F2(k+1:end,1:k,1:k) - permute(reshape(T2(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,1:2),[2,1,3]),[n(2)-k,(n(1)-k)*k]),[2,n(1)-k,k]),[2,1,3]);
uCoeffs(1:k,1:k,1:k) = F1(:,1:k,1:k) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,1:k,1:k),[n(1)-k,k*k]),[2,k,k]);
end
