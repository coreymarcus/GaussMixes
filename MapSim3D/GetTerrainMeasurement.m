function [meas] = GetTerrainMeasurement(m0, bhat, model, terrain)
%GetTerrainMeasurement generates a single terrain measurement

%initialize output
meas = zeros(3,1);

switch model
    case 'Simple'
        
        % Find the range at which the line intersects the terrain
        
        %function for intersect
        fun_inter = @(t) m0(3) + t*bhat(3) - TruthEval(m0(1) + t*bhat(1), m0(2) + t*bhat(2),terrain);
        
        %intersect
        options = optimset('Display','off');
        interrange = fzero(fun_inter,[-10 4*abs(m0(3))],options);
        
        % measurement
        meas = m0 + interrange*bhat;
        
    otherwise
        
        error("Invalid Measurement Model")
        
end


end

