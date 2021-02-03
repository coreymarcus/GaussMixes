function [h] = TruthEval(x, shape)
%TRUTHEVAL Evaluates the true height of the model at a given x value

%locals
n = length(x);
h = zeros(size(x));


% run loop
for ii = 1:n
    
    switch shape
        case 'DoubleRamp'
            
            if(x(ii) < 0)
                h(ii) = x(ii);
            else
                h(ii) = 0.5*x(ii);
            end
        case 'Parabola'
            
            h(ii) = 0.1*x(ii)^2;
            
        otherwise
            disp("Error: Invalid Truth Shape!")
            break
    end
    
end

end

