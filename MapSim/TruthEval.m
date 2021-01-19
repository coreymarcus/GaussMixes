function [h] = TruthEval(x)
%TRUTHEVAL Evaluates the true height of the model at a given x value

%% Simple model
if(x < 0)
    h = -x;
else
    h = 0.5*x;
end

end

