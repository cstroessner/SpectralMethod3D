% Code to generate the plot for Figure 5

addpath('spectral_method_3D')
addpath('tensor_recursive')
clear
clc
close all
rng(1)
tic();
Col = {'b','r','m','g'}; cols = 1;

% initialize params
Ns = 5:2:45; Niters = 1;
for N = Ns
N
iters = 50;
n = [N,N,N];

% intialize function
uZero = @(x,y,z) 1+0.*y.*z.*x;


%% Solve with global spectral method

tic()
% initialize test points
x = rand([100,1])*2-1;
y = rand([100,1])*2-1;
z = rand([100,1])*2-1;

% discretize the functions
u = getFullCoeffsFromFunction(uZero,n);

% discretize the PDE corresponding to L x = u;
% get the Laplacian discretization first
LCP = cell(3);
LCP{1,1} = [0,0,-1];LCP{1,2} = [1,0,0];LCP{1,3} = [1,0,0];
LCP{2,1} = [1,0,0];LCP{2,2} = [0,0,-1];LCP{2,3} = [1,0,0];
LCP{3,1} = [1,0,0];LCP{3,2} = [1,0,0];LCP{3,3} = [0,0,-1];
[opLCP,lambda] = getForwardOperatorUltra(LCP,n);
% add the potential term to opLCP
sinTerm = @(x) sin(pi/2*(x+1));  
V = @(x,y,z) sinTerm(x).*sinTerm(y).*sinTerm(z);
sinCoeffs = get1DCoeffsFromFunction(sinTerm,N);
opLCP{1,4} = MultiplicationMatrix(sinCoeffs,2)*getSUltra(N,1)*getSUltra(N,0);
opLCP{2,4} = MultiplicationMatrix(sinCoeffs,2)*getSUltra(N,1)*getSUltra(N,0);
opLCP{3,4} = MultiplicationMatrix(sinCoeffs,2)*getSUltra(N,1)*getSUltra(N,0);
fastFlag = 0;
% boundary conditions
bc = @(x,y) 0.*x.*y;
[T1,F1,T2,F2,T3,F3] = getBoundaryConditionMatrices(n,bc,bc,bc,bc,bc,bc);

% get mean of the potential
T = getFullCoeffsFromFunction(@(x,y,z) sinTerm(x) .* sinTerm(y) .* sinTerm(z),n);
a = sqrt(L2scalarProduct(T,T));

% get operator for preconditioning
LCPprec = cell(3);
LCPprec{1,1} = [a,0,-1];LCPprec{1,2} = [1,0,0];LCPprec{1,3} = [1,0,0];
LCPprec{2,1} = [1,0,0];LCPprec{2,2} = [0,0,-1];LCPprec{2,3} = [1,0,0];
LCPprec{3,1} = [1,0,0];LCPprec{3,2} = [1,0,0];LCPprec{3,3} = [0,0,-1];
opLCPprec = getForwardOperatorUltra(LCPprec,n);
restart = 10;
maxRestarts = 10;

% inverse iteration
for iter = 1:iters
    
    % step 0 normalize
    u = u./sqrt(L2scalarProduct(u,u));
    uold = u;
    
    % step 1 set up rhs in terms of ultrasphericals
    rhs = tprod(u,getSUltra(n(1),0),getSUltra(n(2),0),getSUltra(n(3),0));
    rhs = tprod(rhs,getSUltra(n(1),1),getSUltra(n(2),1),getSUltra(n(3),1));
    
    % step 2 solve the PDE including the boundary conditions of u (u = solveWithElimination(opLCP,rhs,T1,F1,T2,F2,T3,F3,0);)
    % use preconditioning and GMRES for faste computations
    [opred,rhsred] = getReducedSystem(opLCP,rhs,T1,F1,T2,F2,T3,F3);
    opLCPprecRed = getReducedSystem(opLCPprec,rhs,T1,F1,T2,F2,T3,F3); 
    M = @(x) reshape(solveLinearEquation(opLCPprecRed,reshape(x,n-2),1),[prod(n-2),1]);
    [uCoeffsGMRESredPrec ,~,~,iterGMRESredPrec,resvecGMRESredPrec] = gmres(@(x) reshape( applyForwardOperator(opred,reshape(x,n-2)),[prod(n-2),1]),reshape(rhsred,[prod(n-2),1]),...
        restart,1e-12,maxRestarts,@(x)M(x));
    u = completeReducedSystem(reshape(uCoeffsGMRESredPrec,n-2),2,F1,F2,F3,T1,T2,T3);
    
    % step 3 estimate the eigenvalue
    lambdaGSM(iter) = 1/(L2scalarProduct(uold,u)/L2scalarProduct(uold,uold));
end
tGSM(Niters) = toc();


%% solve with finite differences for comparison
tic()

% discretization params
hx = 2/(N-1);

% discretize initial u
nodes = linspace(-1,1,N);
[X,Y,Z]=ndgrid(nodes,nodes,nodes);
U = uZero(X,Y,Z);
U = U(:);

% assemble differential operator
e = ones(N,1);
% Laplacian
tmpD = (1/hx)^2*spdiags([e -2*e e],-1:1,N,N);
D = sparse(N,N);
D(2:N-1,:) = tmpD(2:N-1,:);
I = speye(N,N);
D = kron(D,kron(I,I)) + kron(I,kron(D,I)) + kron(I,kron(I,D));
D = -D;
% add the potential term
S = spdiags(sinTerm(nodes)',0,N,N);
S = kron(S,kron(S,S));
M = D+S;

% restrict the system to the inner indices to incorporate the homogeneous Dirichlet boundary conditions implicitly 
II = sparse(N,N); II(1,1) = 1; II(end,end) = 1; Inds = kron(II,kron(I,I)) + kron(I,kron(II,I)) + kron(I,kron(I,II)); % now II(i,i) > 0 if boundary point
innerInds = zeros(N^3,1);
for i = 1:N^3
    if Inds(i,i) == 0
        innerInds(i) = 1;
    end
end
innerInds = logical(innerInds);
innerM = M(innerInds,innerInds);

% iterate
for iter = 1:iters
    % step 0 normalize using the L2 norm
    U = U(:);
    U = U./normFD(U,N);
    Uold = U;
    
    % step 1 solve
    U = U(innerInds);
    Utmp = innerM\U;
    
    % step 2 assembe the full U
    U = zeros(N,N,N);
    U = U(:);
    U(innerInds) = Utmp;    
    
    % step 3 estimate the eigenvalues
    lambdaFD(iter) = 1/(prodFD(U,Uold,N)/prodFD(Uold,Uold,N));
end
     
tFD(Niters) = toc();

% store the finial lambdas
lamGSM(Niters) = lambdaGSM(end);
lamFD(Niters) = lambdaFD(end);
Niters = Niters + 1;
end

%% plot convergence rates
close all
set(gca,'fontsize',10)
cols = 1;
set(figure(1), 'Position', [0 0 470 400])
semilogy(Ns(1:end-1),abs(lamFD(1:end-1)-lamFD(end)),'r',Ns(1:end-1),abs(lamGSM(1:end-1)-lamGSM(end)),'b')
xlabel('n','Interpreter','latex')
ylabel('$|\lambda[n] - \lambda[45]|$','Interpreter','latex')
leg = legend('finite differences','global spectral method')
set(leg,'Interpreter','latex','Location','east');
xlim([Ns(1) Ns(end-1)])
%print -depsc 'EigenvalueConvergenceSpeed'

timeGSM = sum(tGSM)
timeFD = sum(tFD)

function val = prodFD(U,V,N)
val = (U'*V)/(N^3)*8;
end
