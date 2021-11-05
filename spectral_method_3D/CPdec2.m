function [U,V,W,cpErr] = CPdec2(T,r)
% Computes the CP decomposition of the tensor T with rank r using the
% Levenberg- Marquardt method.

n = size(T);
U = rand([n(1),r]);
V = rand([n(2),r]);
W = rand([n(3),r]);

[U,V,W] = LevMarMet(T,U,V,W);

TT = CPtoFull(U,V,W);
cpErr = max(abs(T-TT),[],'all')
end

function [A,B,C,iter]=LevMarMet(X,A,B,C)
% Levenberg-Marquardt method for computing CP decompositions
n = size(X); r = size(A,2);
maxIter = 100; iter = 0;
XX = CPtoFull(A,B,C);
Res =X(:)- XX(:);
Jac =Jacobian(A,B,C,n,r);
mu = 0.01*max(diag(Jac'*Jac));
v=2;
a = Jac'*Jac;
g = Jac'*Res;
I = eye(size(a));
tol = 1e-12;
st = [A;B;C];
while norm(Res) >tol && iter<maxIter
    iter = iter + 1;
    h=-pinv(a+mu*I)*g;
    if norm(h)<=tol*(tol+norm(st(:)))
        warning('early break - not enough descent')
        break
    end
    stnew=st(:)+h;
    stnew = reshape(stnew,n(1)+n(2)+n(3),r);
    A = stnew(1:n(1),:);
    B = stnew(n(1)+1:n(1)+n(2),:);
    C = stnew(n(1)+n(2)+1:n(1)+n(2)+n(3),:);
    XX = CPtoFull(A,B,C);
    q=((norm(Res))^2-(norm(XX(:)-X(:)))^2)/(dot(h,mu*h-g));
    if q>0
        st(:)=stnew(:);
        Jac=Jacobian(A,B,C,n,r);
        a=Jac'*Jac;
        Res=X(:)-XX(:);
        g=Jac'*Res;
        mu=mu*max(1/3,1-(2*q-1)^3);
        v=2;
    else
        mu=mu*v;
        v=2*v;
    end
    if iter ==1000
        warning('max iter reached')
        break
    end
end
end

function T = CPtoFull(A,B,C)
r = size(A,2);
n = [size(A,1),size(B,1),size(C,1)];
T = zeros(n);
for i = 1:r
    for a = 1:n(1)
        for b = 1:n(2)
            for c = 1:n(3)
                T(a,b,c) = T(a,b,c) + A(a,i)*B(b,i)*C(c,i);
            end
        end
    end
end
end

function Jac = Jacobian(A,B,C,n,r)
% jacobian for LevMarMet
Jac = zeros([n(1)*n(2)*n(3),(n(1)+n(2)+n(3))*r]);
for I = 1:n(1)*n(2)*n(3)
    [i,j,k] = ind2sub([n(1),n(2),n(3)],I);
    for J = 1:(n(1)+n(2)+n(3))*r
        [a,b] = ind2sub([n(1)+n(2)+n(3),r],J);
        if a == i
            Jac(I,J) = -B(j,b)*C(k,b);
        elseif a-n(1) == j
            Jac(I,J) = -A(i,b)*C(k,b);
        elseif a-n(1)-n(2) == k
            Jac(I,J) = -A(i,b)*B(j,b);
        end
    end
end
end
