function c = evalUltraPolynomials(N,lambda,x)
% returns the function evaluation at x of the first N C^{(lambda)} polynomials
% x can be a row vector of evaluation points 
% c contains the evaluation of x(1) in the first row and x(2) in the second
    c(1,:) = 1*ones(size(x));
    c(2,:) = 2*lambda.*x;
    for i = 1:N-2
        c(i+2,:) = 2*(i+lambda)/(i+1).*x.*c(i+1,:)-(i+2*lambda -1)/(i+1)*c(i,:);
    end
    c = c';
end


