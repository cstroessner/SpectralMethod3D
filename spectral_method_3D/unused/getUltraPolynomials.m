function c = getUltraPolynomials(N,lambda)
% returns function handles to the first N C^{(lambda)} polynomials
% for evaluations evalUltraPolynomials is faster
    c{1} = @(x) 1;
    c{2} = @(x) 2*lambda*x;
    for i = 1:N-2
        c{i+2} = @(x)  2*(i+lambda)/(i+1)*x.*c{i+1}(x)-(i+2*lambda -1)/(i+1)*c{i}(x);
    end
end

