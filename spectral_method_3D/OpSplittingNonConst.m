function [Split,cpErr] = OpSplittingNonConst(L,n,R)
% input: cell array of coefficient functions L{a,b,c} = @(x,y,z)
% alpha(x,y,z), resolution n, desired rank R
% output: cell array of chebyshev coefficients representing the coefficient functions
% of the one-dimensional differential operators in terms of Chebyshev
% coefficients 


N = size(L);

% get the function evalautions for the coefficient functions
for a = 1:3
    for b = 1:3
        for c = 1:3
            alphaEvals{a,b,c} = getEvalsFromFunction(L{a,b,c},n);
        end
    end
end

% construct the big tensor A
A = zeros(N(1),n(1),N(2),n(2),N(3),n(3));
cheb = @(i,n) cos((i-1).*pi/(n-1));
chebs1 = cheb(1:N(1),N(1));
chebs2 = cheb(1:N(2),N(2));
chebs3 = cheb(1:N(3),N(3));
for a = 1:N(1)
    for b = 1:N(2)
        for c = 1:N(3)
            alph = alphaEvals{a,b,c};
            for i = 1:n(1)
                for j = 1:n(2)
                    for k = 1:n(3)
                        for u = 1:N(1)
                            for v = 1:N(2)
                                for w = 1:N(3)
                                    A(u,i,v,j,w,k) = A(u,i,v,j,w,k) + alph(i,j,k) .* chebs1(u).^(a-1).*chebs2(v).^(b-1).*chebs3(w).^(c-1) ;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

T = reshape(A,[n(1)*N(1),n(2)*N(2),n(3)*N(3)]);
%CP decomposition with given rank R
[U,V,W,cpErr] = CPdec2(T,R);


% build output objects
Split = cell(R,3);
for r = 1:R
Ur = U(:,r);
Ur = reshape(Ur,[N(1),n(1)]);
Ur = Vals2ChebCoeffsMat(size(Ur,1))*Ur;
Ur = permute(Vals2ChebCoeffsMat(size(Ur,2))*permute(Ur,[2,1]),[2,1]);
Ur = Ur'*Cheby2Monomial(N(1));

Vr = V(:,r);
Vr = reshape(Vr,[N(2),n(2)]);
Vr = Vals2ChebCoeffsMat(size(Vr,1))*Vr;
Vr = permute(Vals2ChebCoeffsMat(size(Vr,2))*permute(Vr,[2,1]),[2,1]);
Vr = Vr'*Cheby2Monomial(N(2));

Wr = W(:,r);
Wr = reshape(Wr,[N(3),n(3)]);
Wr = Vals2ChebCoeffsMat(size(Wr,1))*Wr;
Wr = permute(Vals2ChebCoeffsMat(size(Wr,2))*permute(Wr,[2,1]),[2,1]);
Wr = Wr'*Cheby2Monomial(N(3));

Split{r,1} = Ur;
Split{r,2} = Vr;
Split{r,3} = Wr;
end

end

function M = Cheby2Monomial(N)
%returns matrix M mapping Chebyshev coefficients to monomial coefficients
M  = zeros(N);
M(1,1) = 1;
if N < 2
    return
end
M(2,2) = 1;
for i = 3:N
M (i,:) = 2*[0,M(i-1,1:end-1)] -M(i-2,:);
end
end


