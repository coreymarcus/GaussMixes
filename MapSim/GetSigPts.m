function [Xi,w] = GetSigPts(mu, P, k)
%GetSigPts returns the points, Xi, and Weights, w, extracted from a gaussian
%according to the 2n+1 procedure used for the UKF

%length of state
n = length(mu);

%initialize
Xi = zeros(n,2*n+1);
w = zeros(1,2*n+1);

%Matrix square root
S = chol(P,'lower');

for ii = 1:n
    Xi(:,ii) = mu - sqrt(n+k)*S(:,ii);
    Xi(:,ii + n) = mu + sqrt(n+k)*S(:,ii);
    w(ii) = .5/(n+k);
    w(ii+n) = .5/(n+k);
end

%center point
Xi(:,2*n+1) = mu;
w(2*n+1) = k/(n+k);

end

