function [uCoeffs,ResidualNorm] = solvePDE(LCP,rhsFunction,rightbc,leftbc,topbc,bottombc,frontbc,backbc)
% a function wrapper which solves a PDE

% compute the CP decomposition if the coefficient tensor is given
if ndims(LCP) == 3 
    LCP = CPdec(LCP);
end
% set the fast Flag if fast solving is possible
[LCP,fastFlag] = checkStructure(LCP);

% TODO fix n and refine
n = [6, 6, 6]; ResidualNorm = Inf;
while ResidualNorm > 1e-10 %better check if the individual fibers are converged
n(:) = floor(sqrt(2).^(floor(2*log2(n(:))) + 1)) + 1; % same as in Chebfun3 Phase 

% solves the pde
[opLCP,lambda] = getForwardOperatorUltra(LCP,n); 
rhsCoeffs = getFullCoeffsFromFunctionUltra(rhsFunction,n,lambda);

% up to here optimized for sparsity and arbitrary ranks

[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc);
uCoeffs = solveWithElimination(opLCP,rhsCoeffs,T1,F1,T2,F2,T3,F3,fastFlag);

R1 = max(abs(tprod(uCoeffs,T1,eye(n(2)),eye(n(3)))-F1),[],'all');
R2 = max(abs(tprod(uCoeffs,eye(n(1)),T2,eye(n(3)))-F2),[],'all');
R3 = max(abs(tprod(uCoeffs,eye(n(1)),eye(n(2)),T3)-F3),[],'all');
Rop =  max(abs(applyForwardOperator(opLCP,uCoeffs) - rhsCoeffs),[],'all');
ResidualNorm = max([R1,R2,R3,Rop]);

if (n(1) > 17 && fastFlag == 0) || (n(1) > 1250)
    warning('Maximum for grid resolution n reached!')
    return
end
end