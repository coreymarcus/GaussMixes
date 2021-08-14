function [phi] = UnitVectorLog(n)
%UnitVectorLog 

%locals
kappa = n(1:2);
lambda = n(3);
kappa_norm = norm(kappa);
kappa_hat = kappa/kappa_norm;

%check length
if(abs(norm(n) - 1) > 0.0001)
    warning('Non-unit vector input')
end

if(kappa_norm == 0)
    phi = zeros(2,1);
else
    phi = atan2(kappa_norm,lambda)*kappa_hat;
end

end

