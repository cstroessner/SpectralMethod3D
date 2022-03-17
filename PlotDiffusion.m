% code to generate the plot in Figure 3
addpath('spectral_method_3D')
addpath('tensor_recursive')
clear
close all
rng(1)

% inizialize parameters
a1 = @(x) x.^2+1; a2 = @(y) y.^2+1; a3 = @(z) z.^2+1;
b1 = @(x) exp(x); b2 = @(y) exp(y); b3 = @(z) exp(z);
a = @(x,y,z) a1(x).*a2(y).*a3(z)+b1(x).*b2(y).*b3(z);
utrue = @(x,y,z)  sin(pi*x).*sin(pi*y).*sin(pi*z);
f = @(x,y,z) -( ...
    pi*sin(pi*y).*sin(pi.*z) .* (2*x.*cos(pi*x) - pi*a1(x).*sin(pi*x)).*a2(y).*a3(z) + ...
    pi*sin(pi*x).*sin(pi.*z) .* (2*y.*cos(pi*y) - pi*a2(y).*sin(pi*y)).*a1(x).*a3(z) + ...
    pi*sin(pi*y).*sin(pi.*x) .* (2*z.*cos(pi*z) - pi*a3(z).*sin(pi*z)).*a1(x).*a2(y) + ...
    pi*sin(pi*y).*sin(pi.*z) .* (exp(x).*cos(pi*x) - pi*b1(x).*sin(pi*x)).*b2(y).*b3(z) + ...
    pi*sin(pi*x).*sin(pi.*z) .* (exp(y).*cos(pi*y) - pi*b2(y).*sin(pi*y)).*b1(x).*b3(z) + ...
    pi*sin(pi*y).*sin(pi.*x) .* (exp(z).*cos(pi*z) - pi*b3(z).*sin(pi*z)).*b1(x).*b2(y) );

% discretization
N = 31;
n = [N,N,N];
% operator
opLCPa = getForwardOperatorDiffusionProblem(a1,a2,a3,n);
[opLCPb,lambda] = getForwardOperatorDiffusionProblem(b1,b2,b3,n);
for mode = 1:3
    for i = 1:3
        opLCP{mode,i} = opLCPa{mode,i};
    end
    for i = 4:6
        opLCP{mode,i} = opLCPb{mode,i-3}; 
    end
end
fastFlag = 0;
% righ hand side
rhsCoeffs = getFullCoeffsFromFunction(f,n);
rhsCoeffsUltra = getFullCoeffsFromFunctionUltra(f,n,lambda);
% boundary conditions
rightbc = @(y,z) utrue(1,y,z);
leftbc = @(y,z) utrue(-1,y,z);
topbc = @(x,z) utrue(x,1,z);
bottombc = @(x,z) utrue(x,-1,z);
frontbc = @(x,y) utrue(x,y,-1);
backbc = @(x,y) utrue(x,y,1);
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,rightbc,leftbc,topbc,bottombc,frontbc,backbc);

% build preconditioner (constant)
T = getFullCoeffsFromFunction(a,n);
const = sqrt(L2scalarProduct(T,T));
c1 = @(x) 0.*x+const; c2 = @(y) 0.*y+1; c3 = @(z) 0.*z+1;
opLCPprec = getForwardOperatorDiffusionProblem(c1,c2,c3,n);
opLCPprecRed = getReducedSystem(opLCPprec,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);

% build preconditioner (rank-1)
opLCPprec2 = getForwardOperatorDiffusionProblem(a1,a2,a3,n);
opLCPprecRed2 = getReducedSystem(opLCPprec2,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);
fastFlagPrec = 1; restart = 15; maxRestarts = 15;

% solve via reshaping
tic()
uReshape = solveWithElimination(opLCP,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3,0);
tReshape = toc()

% GMRES no preconditioning
tic()
[opred,rhsred] = getReducedSystem(opLCP,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);
[ured ,~,~,itersNoPrec,resNoPrec] = gmres(@(x) reshape( applyForwardOperator(opred,reshape(x,n-2)),[prod(n-2),1]),reshape(rhsred,[prod(n-2),1]),...
    restart,1e-12,maxRestarts);
uGMRESnoPrec = completeReducedSystem(reshape(ured,n-2),2,F1,F2,F3,T1,T2,T3);
tGMRESnoPrec = toc();

% GMRES with preconditioning (constant)
tic()
[opred,rhsred] = getReducedSystem(opLCP,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);
M = @(x) reshape(solveLinearEquation(opLCPprecRed,reshape(x,n-2),1),[prod(n-2),1]);
[ured ,~,~,itersWithPrec,resWithPrec] = gmres(@(x) reshape( applyForwardOperator(opred,reshape(x,n-2)),[prod(n-2),1]),reshape(rhsred,[prod(n-2),1]),...
    restart,1e-12,maxRestarts,@(x)M(x));
uGMRESwithPrec = completeReducedSystem(reshape(ured,n-2),2,F1,F2,F3,T1,T2,T3);
tGMRESwithPrec = toc()

% GMRES with preconditioning (rank-1)
tic()
[opred,rhsred] = getReducedSystem(opLCP,rhsCoeffsUltra,T1,F1,T2,F2,T3,F3);
M = @(x) reshape(solveLinearEquation(opLCPprecRed2,reshape(x,n-2),1),[prod(n-2),1]);
[uredR1 ,~,~,itersWithPrecR1,resWithPrecR1] = gmres(@(x) reshape( applyForwardOperator(opred,reshape(x,n-2)),[prod(n-2),1]),reshape(rhsred,[prod(n-2),1]),...
    restart,1e-12,maxRestarts,@(x)M(x));
uGMRESwithPrecR1 = completeReducedSystem(reshape(uredR1,n-2),2,F1,F2,F3,T1,T2,T3);
tGMRESwithPrecR1 = toc()

% error analysis
uTrue = getFullCoeffsFromFunction(utrue,n);
errCoeffsReshape = max(abs(uReshape-uTrue),[],'all');
errCoeffsNoPrec = max(abs(uGMRESnoPrec-uTrue),[],'all');
errCoeffsWithPrec = max(abs(uGMRESwithPrec-uTrue),[],'all');
errCoeffsWithPrecR1 = max(abs(uGMRESwithPrecR1-uTrue),[],'all');
for iters = 1:1000
    x = 2*(rand(1)-0.5); y = 2*(rand(1)-0.5); z = 2*(rand(1)-0.5);
    errReshape(iters) = abs(funeval(uReshape,x,y,z)-utrue(x,y,z));
    errNoPrec(iters) = abs(funeval(uGMRESnoPrec,x,y,z)-utrue(x,y,z));
    errWithPrec(iters) = abs(funeval(uGMRESwithPrec,x,y,z)-utrue(x,y,z));
    errWithPrecR1(iters) = abs(funeval(uGMRESwithPrecR1,x,y,z)-utrue(x,y,z));
    errInterpolation(iters) = abs(funeval(uTrue,x,y,z)-utrue(x,y,z));
end
errReshape = max(errReshape);
errNoPrec = max(errNoPrec);
errWithPrec = max(errWithPrec);
errWithPrecR1 = max(errWithPrecR1);
errInterpolation = max(errInterpolation);

%% plot 
set(gca,'fontsize',10)
cols = 1;
set(figure(1), 'Position', [0 0 470 400])
semilogy(1:numel(resNoPrec),resNoPrec,'r',1:numel(resWithPrec),resWithPrec,'b',1:numel(resWithPrecR1),resWithPrecR1,'m')
xlabel('GMRES iteration')
ylabel('residual')
cols = cols+1;
leg = legend('no preconditioning','$b(x,y,z) = ||a||_{\mathcal{L}^2}$','$b(x,y,z) = (1+x^2)(1+y^2)(1+z^2)$');
set(leg,'Interpreter','latex','Location','northeast');
xlim([1,100])
print -depsc 'GMRESresidualVsIterations'