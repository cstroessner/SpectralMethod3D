% code to reproduce the plots in Section 5.2.5
addpath('spectral_method_3D')
addpath('tensor_recursive')
addpath('tensorlab_2016-03-28')
clear
close all
rng(1)

% solution and coefficient
k = @(x,y,z) sqrt(x+y+z+42);

%discretization
sizes = 45:-5:5;
for ITER = 1:numel(sizes)
N = sizes(ITER)+1
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
f = @(x,y,z) 1+0.*x.*y.*z;
% boundary conditions
rightbc = @(y,z) 0.*y.*z;
leftbc = @(y,z) 0.*y.*z;
topbc = @(x,z) 0.*x.*z;
bottombc = @(x,z) 0.*x.*z;
frontbc = @(x,y) 0.*x.*y;
backbc = @(x,y) 0.*x.*y;

% CP decomposition and operator splitting
tic();
[Split,cpErr] = OpSplittingNonConst(L,n,3+8);
tCPdec(ITER) = toc();

% discretize
[opLCP,lambda] = getForwardOperatorUltraNonConst(Split,n);
rhsCoeffsUltra = getFullCoeffsFromFunctionUltra(f,n,lambda);
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc);
% change Dirichlet to Neumann boundary conditions on the right 
T1 = zeros(size(T1));
for i = 1:size(T1,2)
   T1(1,i) = (-1)^(i-1);
   T1(2,i) = (i-1)^2;
end
T1 = T1(1:2,1:2)\T1;

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
tGMRESwithPrec(ITER) = toc();

if ITER == 1
   uControl = uGMRESwithPrec;
end

for iters = 1:1000
    x = 2*(rand(1)-0.5); y = 2*(rand(1)-0.5); z = 2*(rand(1)-0.5);
    uErrGMRES(iters) = abs(funeval(uGMRESwithPrec,x,y,z)-funeval(uControl,x,y,z));
end

estErr(ITER) = max(uErrGMRES);
CPerrs(ITER) = cpErr;
residuals(ITER) = computeResidual(uGMRESwithPrec,opLCP,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);
ITER = ITER +1;
end

%% plot
close all
set(gca,'fontsize',10)
set(figure(1), 'Position', [0 0 470 400])
semilogy(sizes(2:end),residuals(2:end),'k--')
hold on
semilogy(sizes(2:end),CPerrs(2:end),'k:') 
semilogy(sizes(2:end),estErr(2:end),'k-') 
xlabel('n')
leg = legend('residual','CP-approx. error','$||u[n]-u[45]||_\infty$');
set(leg,'Interpreter','latex','Location','northeast');
print -depsc 'HelmholtzUnknownDecays'

%% plot 2
figure(2)
[X,Y] = meshgrid(-1:0.05:1,-1:0.05:1);
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Z(i,j) = funeval(uControl,X(i,j),Y(i,j),1/4);
    end
end
s = surf(X,Y,Z,'FaceAlpha',0.5)
v = [-5 -3 5];
[caz,cel] = view(v)
%print -depsc 'HelmholtzSolutionUnknown'