function [meas] = GetTerrainMeasurement(x, theta, R, model, terrain)
%GetTerrainMeasurement generates a single terrain measurement given pose x,
%measurement bearing theta, covariance R, and model

%initialize output
meas = zeros(2,1);

switch model  
    case 'Simple'
        
        %angle
        alpha = x(3) + theta;
        
        %find equation for line at angle alpha and where it intersects
        %terrain
        if(rem(alpha,pi/2) == 0)
            
            %measurement is straight down so this is easy
            meas(1) = x(1);
            meas(2) = TruthEval(x(1),terrain);
            
        else
            
            %now measurement is along a line so we should solve for the
            %intersection
            
            %slope
            m = tan(alpha);
            
            %intercept
            b = x(2) - m*x(1);
            
            %function for intersect
            fun_inter = @(x_eval) m*x_eval + b - TruthEval(x_eval,terrain);
            
            if(m < 0)
                bounds = [x(1) 1000];
            else
                bounds = [-1000 x(1)];
            end
            
            %intersect
            meas(1) = fzero(fun_inter,bounds);
            meas(2) = TruthEval(meas(1),terrain);    
            
        end
        
        
    otherwise
        
        disp("Error: Invalid Measurement Model")
        
end
        
        
end

