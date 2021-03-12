function [slope] = GetGMSlopeEstimate(gauss_list,x_eval, method)
%GetGMSlopeEstimate does what it says

%locals
Ngauss = length(gauss_list);
slopeweight = zeros(Ngauss,1);
slopes = zeros(Ngauss,1);

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
    slopeweight(ii) = obj.GaussEvalST(x_eval, y_eval) / obj.P_mb(1,1);
    
    %storage for later
    slopes(ii) = m_hat;
    
end

%get slope
switch method
    case "ML"
        [~, maxidx] = max(slopeweight);
        slope = slopes(maxidx(1));
        
    case "MMSE"
        
        %normalize weights
        slopeweight = slopeweight./sum(slopeweight);
        slope = sum(slopeweight.*slopes);
        
    otherwise
        disp("Error: Invalid Slope Method!")
end


end

