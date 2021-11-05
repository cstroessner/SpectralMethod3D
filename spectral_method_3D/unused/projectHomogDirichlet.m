% f = @(x,y,z) (x-1).*(x+1).*(y-1).*(y+1).*(z-1).*(z+1);
% uCoeffsUlt = getFullCoeffsFromFunctionUltra(f,[10,10,10],[2,2,2]);
% uCoeffs = getFullCoeffsFromFunction(f,[10,10,10]);
% uCoeffs2 = projectHomogDirichlet1(uCoeffsUlt);
% 
% err = max(abs(uCoeffs2-uCoeffs),[],'all')
%sanityCheck = max(abs(tprod(uCoeffsUlt,eye(10),evalUltraPolynomials(10,2,[1,-1]),eye(10))),[],'all')


function uCoeffs = projectHomogDirichlet(uCoeffs)
% projects a function given in terms of C 2 ultraspherical coefficients onto a
% function with zero boundary conditions
% afterwards the coefficients are mapped to chebyshev coefficients

lambda = 2; k = 2;
n = size(uCoeffs);
C1 = evalUltraPolynomials(n(1),lambda,[1,-1]);
Cprincsub1 = C1(end-1:end,end-1:end);
C1 = Cprincsub1\C1;
C2 = evalUltraPolynomials(n(2),lambda,[1,-1]);
Cprincsub2 = C2(end-1:end,end-1:end);
C2 = Cprincsub2\C2;
C3 = evalUltraPolynomials(n(3),lambda,[1,-1]);
Cprincsub3 = C3(end-1:end,end-1:end);
C3 = Cprincsub3\C3;
T1 = C1;
T2 = C2;
T3 = C3;

% no change of F terms since they are zero
F1 = zeros([2,n(2),n(3)]);
F2 = zeros([n(1),2,n(3)]);
F3 = zeros([n(1),n(2),2]);
F1 = reshape( Cprincsub1\reshape(F1,[2,n(2)*n(3)]),[2,n(2),n(3)]);
F2 = permute(reshape(Cprincsub2\reshape(permute(F2,[2,1,3]),[2,n(1)*n(3)]),[2,n(1),n(3)]),[2,1,3]);
F3 = permute(reshape(Cprincsub3\reshape(permute(F3,[3,2,1]),[2,n(2)*n(1)]),[2,n(2),n(1)]),[3,2,1]);

% this changes the first 2 coefficients - I should change the last 2 coefficients in each mode
% uCoeffs(1:2,k+1:end,k+1:end) = F1(:,k+1:end,k+1:end) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,k+1:end,k+1:end),[n(1)-k,(n(2)-k)*(n(3)-k)]),[2,n(2)-k,n(3)-k]);
% uCoeffs(k+1:end,1:2,k+1:end) = F2(k+1:end,:,k+1:end) - permute(reshape(T2(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,k+1:end),[2,1,3]),[n(2)-k,(n(1)-k)*(n(3)-k)]),[2,n(1)-k,n(3)-k]),[2,1,3]);
% uCoeffs(k+1:end,k+1:end,1:2) = F3(k+1:end,k+1:end,:) - permute(reshape(T3(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,k+1:end),[3,2,1]),[n(3)-k,(n(2)-k)*(n(1)-k)]),[2,n(2)-k,n(1)-k]),[3,2,1]);
% uCoeffs(1:k,k+1:end,1:k) = F1(:,k+1:end,1:k) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,k+1:end,1:2),[n(1)-k,(n(2)-k)*k]),[2,n(2)-k,k]);
% uCoeffs(1:k,1:k,k+1:end) = F1(:,1:k,k+1:end) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,1:2,k+1:end),[n(1)-k,(n(3)-k)*k]),[2,k,n(3)-k]);
% uCoeffs(k+1:end,1:k,1:k) = F2(k+1:end,1:k,1:k) - permute(reshape(T2(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,1:2),[2,1,3]),[n(2)-k,(n(1)-k)*k]),[2,n(1)-k,k]),[2,1,3]);
% uCoeffs(1:k,1:k,1:k) = F1(:,1:k,1:k) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,1:k,1:k),[n(1)-k,k*k]),[2,k,k]);

% new approach to change the last 2 coefficients in each mode
uCoeffs(end-1:end,1:end-2,1:end-2) = F1(:,1:end-2,1:end-2) - reshape(T1(:,1:end-2)*reshape(uCoeffs(1:end-2,1:end-2,1:end-2),[n(1)-k,(n(2)-k)*(n(3)-k)]),[2,n(2)-k,n(3)-k]);
uCoeffs(:,end-1:end,1:end-2) = F2(:,:,1:end-2) - permute(reshape(T2(:,1:end-2)*reshape(permute(uCoeffs(:,1:end-2,1:end-2),[2,1,3]),[n(2)-k,(n(1))*(n(3)-k)]),[2,n(1),n(3)-k]),[2,1,3]);
uCoeffs(:,:,end-1:end) = F3(:,:,:) - permute(reshape(T3(:,1:end-2)*reshape(permute(uCoeffs(:,:,1:end-2),[3,2,1]),[n(3)-k,(n(2))*(n(1))]),[2,n(2),n(1)]),[3,2,1]);

% C = evalUltraPolynomials(n(1),lambda,[1,-1]);
% err = max( [max(abs(tprod(uCoeffs,eye(10),C,eye(10))-F2),[],'all'), ...
%     max(abs(tprod(uCoeffs,eye(10),eye(10),C)-F3),[],'all'), ...
%     max(abs(tprod(uCoeffs,C,eye(10),eye(10))-F1),[],'all')])

% map to cheby coeffs 
uCoeffs = invtprod(uCoeffs,getSUltra(n(1),1),getSUltra(n(2),1),getSUltra(n(3),1));
uCoeffs = invtprod(uCoeffs,getSUltra(n(1),0),getSUltra(n(2),0),getSUltra(n(3),0));
end

