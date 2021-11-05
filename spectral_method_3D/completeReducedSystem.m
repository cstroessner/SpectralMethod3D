function uCoeffs = completeReducedSystem(uCoeffs,k,F1,F2,F3,T1,T2,T3)
% Given U_222 as solution of the reduced system obtained after
% substitution, this function adds the the remaining blocks to reconstruct
% U form the boundary conditions

n = size(uCoeffs)+k;
uCoeffs2 = zeros(n);
uCoeffs2(k+1:end,k+1:end,k+1:end) = uCoeffs;
uCoeffs = uCoeffs2;

% find the other blocks by using the boundary conditions
uCoeffs(1:2,k+1:end,k+1:end) = F1(:,k+1:end,k+1:end) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,k+1:end,k+1:end),[n(1)-k,(n(2)-k)*(n(3)-k)]),[2,n(2)-k,n(3)-k]);
uCoeffs(k+1:end,1:2,k+1:end) = F2(k+1:end,:,k+1:end) - permute(reshape(T2(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,k+1:end),[2,1,3]),[n(2)-k,(n(1)-k)*(n(3)-k)]),[2,n(1)-k,n(3)-k]),[2,1,3]);
uCoeffs(k+1:end,k+1:end,1:2) = F3(k+1:end,k+1:end,:) - permute(reshape(T3(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,k+1:end),[3,2,1]),[n(3)-k,(n(2)-k)*(n(1)-k)]),[2,n(2)-k,n(1)-k]),[3,2,1]);
uCoeffs(1:k,k+1:end,1:k) = F1(:,k+1:end,1:k) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,k+1:end,1:2),[n(1)-k,(n(2)-k)*k]),[2,n(2)-k,k]);
uCoeffs(1:k,1:k,k+1:end) = F1(:,1:k,k+1:end) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,1:2,k+1:end),[n(1)-k,(n(3)-k)*k]),[2,k,n(3)-k]);
uCoeffs(k+1:end,1:k,1:k) = F2(k+1:end,1:k,1:k) - permute(reshape(T2(:,k+1:end)*reshape(permute(uCoeffs(k+1:end,k+1:end,1:2),[2,1,3]),[n(2)-k,(n(1)-k)*k]),[2,n(1)-k,k]),[2,1,3]);
uCoeffs(1:k,1:k,1:k) = F1(:,1:k,1:k) - reshape(T1(:,k+1:end)*reshape(uCoeffs(k+1:end,1:k,1:k),[n(1)-k,k*k]),[2,k,k]);
end