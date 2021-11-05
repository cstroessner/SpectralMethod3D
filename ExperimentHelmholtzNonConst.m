% code to reproduce the numbers given in Section 5.2.4
addpath('spectral_method_3D')
addpath('tensor_recursive')
clear
close all
rng(1)

% solution and coefficient
utrue = @(x,y,z)  sin(pi*x).*sin(pi*y).*sin(pi*z);
k = @(x,y,z) cos(x+y+z);

%discretization
N = 30;
n = [N,N,N];
%operator
for a = 1:3
    for b = 1:3
        for c = 1:3
            L{a,b,c} = @(x,y,z) 0.*x.*y.*z;
        end
    end
end
L{3,1,1} = @(x,y,z) 1+0.*x.*y.*z;
L{1,3,1} = @(x,y,z) 1+0.*x.*y.*z;
L{1,1,3} = @(x,y,z) 1+0.*x.*y.*z;
L{1,1,1} = @(x,y,z) k(x,y,z);
% rhs
f = @(x,y,z) -3*pi*pi*utrue(x,y,z) + k(x,y,z).*utrue(x,y,z);
% boundary conditions
rightbc = @(y,z) utrue(1,y,z);
leftbc = @(y,z) utrue(-1,y,z);
topbc = @(x,z) utrue(x,1,z);
bottombc = @(x,z) utrue(x,-1,z);
frontbc = @(x,y) utrue(x,y,-1);
backbc = @(x,y) utrue(x,y,1);

% CP decomposition and operator splitting
tic();
Split = OpSplittingNonConst(L,n,6);
tCPdec = toc()

% discretize
[opLCP,lambda] = getForwardOperatorUltraNonConst(Split,n);
rhsCoeffsUltra = getFullCoeffsFromFunctionUltra(f,n,lambda);
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc);

% build preconditioner  
LCPprec = cell(3);
LCPprec{1,1} = [0,0,1];
LCPprec{1,2} = [1,0,0];
LCPprec{1,3} = [1,0,0];
LCPprec{2,1} = [1,0,0];
LCPprec{2,2} = [0,0,1];
LCPprec{2,3} = [1,0,0];
LCPprec{3,1} = [1,0,0];
LCPprec{3,2} = [1,0,0];
LCPprec{3,3} = [0,0,1];
opLCPprec = getForwardOperatorUltra(LCPprec,n);
T = getFullCoeffsFromFunction(k,n);
const = sqrt(L2scalarProduct(T,T));
opLCPprec{1,1} = opLCPprec{1,1} + const*getSUltra(n(1),1)*getSUltra(n(1),0);
opLCPprecRed = getReducedSystem(opLCPprec,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);
fastFlagPrec = 1; restart = 15; maxRestarts = 15;

% GMRES with preconditioning (constant)
tic()
[opred,rhsred] = getReducedSystem(opLCP,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);
M = @(x) reshape(solveLinearEquation(opLCPprecRed,reshape(x,n-2),1),[prod(n-2),1]);
[ured ,~,~,itersWithPrec,resWithPrec] = gmres(@(x) reshape( applyForwardOperator(opred,reshape(x,n-2)),[prod(n-2),1]),reshape(rhsred,[prod(n-2),1]),...
    restart,1e-12,maxRestarts,@(x)M(x));
uGMRESwithPrec = completeReducedSystem(reshape(ured,n-2),2,F1,F2,F3,T1,T2,T3);
tGMRESwithPrec = toc()

% solve with reshape
tic()
uReshape = solveWithElimination(opLCP,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3,0);
tReshape = toc()

% error analysis
uInterp = getFullCoeffsFromFunction(utrue,n);
for iters = 1:1000
    x = 2*(rand(1)-0.5); y = 2*(rand(1)-0.5); z = 2*(rand(1)-0.5);
    uErrReshape(iters) = abs(funeval(uReshape,x,y,z)-utrue(x,y,z));
    uErrGMRES(iters) = abs(funeval(uGMRESwithPrec,x,y,z)-utrue(x,y,z));
    interpErr(iters) = abs(funeval(uInterp,x,y,z)-utrue(x,y,z));
end
uErrReshape = max(uErrReshape)
uErrGMRES = max(uErrGMRES)
interpErr = max(interpErr)

% time small CP decomposition for comparison
T = getFullCoeffsFromFunction(k,n);
tic()
CPdec2(T,3);
tCPsmall = toc()