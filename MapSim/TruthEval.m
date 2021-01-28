function [h] = TruthEval(x)
%TRUTHEVAL Evaluates the true height of the model at a given x value

%locals
n = length(x);
h = zeros(size(x));


% run loop
for ii = 1:n
    
    %% Simple model
    
    if(x(ii) < 0)
        h(ii) = x(ii);
    else
        h(ii) = 0.5*x(ii);
    end
    
end

end

