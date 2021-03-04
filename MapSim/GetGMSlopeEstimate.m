function [slope] = GetGMSlopeEstimate(gauss_list,x_eval)
%GetGMSlopeEstimate does what it says

%locals
Ngauss = length(gauss_list);
slopeweight = zeros(Ngauss,1);

%loop
for ii = 1:Ngauss
    
    %target
    obj = gauss_list{ii};
    
    %estimates
    m_hat = obj.mu_mb(1);
    b_hat = obj.mu_mb(2);
    
    %y coord at this x-value
    y_eval = m_hat*x_eval + b_hat;
    
    %weight
    slopeweight(ii) = obj.Nobs * obj.GaussEvalST(x_eval, y_eval) / obj.P_mb(1,1);
    
end

%get slope
[~, maxidx] = max(slopeweight);
slope = gauss_list{maxidx(1)}.mu_mb(1);

end

