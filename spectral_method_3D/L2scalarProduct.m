function val = L2scalarProduct(U,V)
% this function takes 2 chebyshev coefficient tensors of the same size and
% returns an approximation of the L2 scalar product of the two
% corresponding functions
n = size(U);
F1 = Vals2ChebCoeffsMat(n(1));
F2 = Vals2ChebCoeffsMat(n(2));
F3 = Vals2ChebCoeffsMat(n(3));
UV = tprod(invtprod(U,F1,F2,F3).*invtprod(V,F1,F2,F3),F1,F2,F3);
IntTs1  = getIntegralValues(n(1));
IntTs2  = getIntegralValues(n(2));
IntTs3  = getIntegralValues(n(3));
val = UV;
val = IntTs1'*reshape(val,n(1),n(2)*n(3));
val = IntTs2'*reshape(val,n(2),n(3))*IntTs3;
end

function IntTs  = getIntegralValues(n) 
% see Theorem 19.2 in Trefethens book on approximation theory
    IntTs = zeros(n,1);
    for k = 0:n-1
        if mod(k,2) == 0
            IntTs(k+1) = 2/(1-k^2);
        end        
    end
end
