function [h, slope] = TruthEval(x, shape)
%TRUTHEVAL Evaluates the true height of the model at a given x value

%locals
n = length(x);
h = zeros(size(x));
slope = zeros(size(x));

% run loop
for ii = 1:n
    
    switch shape
        case 'DoubleRamp'
            
            if(x(ii) < 0)
                h(ii) = x(ii);
                slope(ii) = 1;
            else
                h(ii) = 0.5*x(ii);
                slope(ii) = 0.5;
            end
        case 'Parabola'
            
            h(ii) = 0.05*x(ii)^2;
            slope(ii) = 0.10*x(ii);
            
        case 'Rock'
            if(abs(x(ii)) > 0.5)
                h(ii) = 0;
                slope(ii) = 0;
            else
                h(ii) = cos(pi*x(ii));
                slope(ii) = -pi*sin(pi*x(ii));
            end
            
        case 'Line'
            m = 1;
            b = 2;
            h(ii) = m*x(ii) + b;
            slope(ii) = m;
            
        otherwise
            disp("Error: Invalid Truth Shape!")
            break
    end
    
end

end

