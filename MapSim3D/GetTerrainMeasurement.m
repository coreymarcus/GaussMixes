function [meas] = GetTerrainMeasurement(m0, bhat, model, terrain)
%GetTerrainMeasurement generates a single terrain measurement given pose x,
%measurement bearing theta, covariance R, and model

%initialize output
meas = zeros(3,1);

switch model
    case 'Simple'
        
        % Find the range at which the line intersects the terrain
        
        %function for intersect
        fun_inter = @(t) m0(3) + t*bhat(3) - TruthEval(m0(1) + t*bhat(1), m0(2) + t*bhat(2),terrain);
        
        %intersect
        interrange = fzero(fun_inter,[0 10000]);
        
        % measurement
        meas = m0 + interrange*bhat;
        
    otherwise
        
        disp("Error: Invalid Measurement Model")
        
end


end

