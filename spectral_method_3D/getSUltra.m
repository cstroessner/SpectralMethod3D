function S = getSUltra(n,lambda)
% Maps C^{lambda} Coefficients to C^{lambda+1} Coefficients. See S_{lambda} 
% in Townsend and Olver 2015.
if lambda == 0
    v0 = 0.5*ones([n,1]);
    v0(1) = 1;
    v2 = -0.5*ones([n-2,1]);
    v2ext = zeros([n,1]);
    v2ext(3:end) = v2;
    V = [v0,v2ext];
    S = spdiags(V,[0,2],n,n);
    return
end
v = lambda:1:lambda+n-1;
v0 = lambda./v;
v2 = -lambda./v(3:end);
v2ext = zeros([n,1]);
v2ext(3:end) = v2;
V = [v0',v2ext];
S = spdiags(V,[0,2],n,n);
end