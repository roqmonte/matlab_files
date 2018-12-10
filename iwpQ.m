function out = iwpQ(v,ixpx)
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: draws the covariance matrix from the inverse Wishart distribution
% for scale parameter v.
k  = size(ixpx,1);
z  = zeros(v,k);
mu = zeros(k,1);

for i = 1:v
    z(i,:) = (chol(ixpx)'*randn(k,1))';
end

out=inv(z'*z);

