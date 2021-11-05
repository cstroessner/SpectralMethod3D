function [opLCP,D] = getReducedSystem(opLCP,D,T1,F1,T2,F2,T3,F3)
% This function sets up all terms for the reduced system obtained via
% substitution

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
opLCP{1,i} = opLCP{1,i} - opLCP{1,i}(:,1:k)*T1;
opLCP{1,i} = opLCP{1,i}(1:end-k,k+1:end);
opLCP{2,i} = opLCP{2,i} - opLCP{2,i}(:,1:k)*T2;
opLCP{2,i} = opLCP{2,i}(1:end-k,k+1:end);
opLCP{3,i} = opLCP{3,i} - opLCP{3,i}(:,1:k)*T3;
opLCP{3,i} = opLCP{3,i}(1:end-k,k+1:end);
end