function res = computeResidual(uCoeffs,opLCP,rhs,T1,F1,T2,F2,T3,F3)
% computes the residual

n = size(uCoeffs);
R1 = max(abs(tprod(uCoeffs,T1,eye(n(2)),eye(n(3)))-F1),[],'all');
R2 = max(abs(tprod(uCoeffs,eye(n(1)),T2,eye(n(3)))-F2),[],'all');
R3 = max(abs(tprod(uCoeffs,eye(n(1)),eye(n(2)),T3)-F3),[],'all');
Rop =  max(abs(applyForwardOperator(opLCP,uCoeffs) - rhs),[],'all');
res = max([R1,R2,R3,Rop]);
end