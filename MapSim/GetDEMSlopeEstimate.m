function [slope] = GetDEMSlopeEstimate(xhatDEM, grid, x_eval)
%GetDEMSlopeEstimate does what it says

%sort the eval point
idx = discretize(x_eval,grid);

%assume constant size grid
dx = grid(2) - grid(1);

if(idx == 1) %use only forward point
    
    y0 = xhatDEM(idx);
    yR = xhatDEM(idx+1);
    slope = (yR - y0)/dx;
    
elseif(idx == (length(grid) - 1)) %use only backward point
    
    y0 = xhatDEM(idx);
    yL = xhatDEM(idx-1);
    slope = (y0 - yL)/dx;
    
else %average forward and backward
    
    y0 = xhatDEM(idx);
    yR = xhatDEM(idx+1);
    yL = xhatDEM(idx-1);
    slope = 0.5*(yR - y0)/dx + 0.5*(y0 - yL)/dx;
    
    
    
end


end

