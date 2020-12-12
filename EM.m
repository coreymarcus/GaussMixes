function [EMmodel] = EM(model, y, options)
%EM performs some simple expectation maximization

%extract local variables
k = model.k; %number of gaussians
mu_hat = model.mu; %means
var_hat = model.var; %variances
w_hat = model.w; %weights
N = length(y); %number of measurements
gamma = zeros(N,k); %prob that (row) measurement from (col) gaussian

% EM parameters
loopmax = options.loopmax;
tol = options.tol;

% compute the initial log liklihood
L = 0;
for ii = 1:N
    Liter = 0;
    for jj = 1:k
        Liter = Liter + w_hat(jj)*gaussEval(y(ii),mu_hat(jj),var_hat(jj));
    end
    L = L + log(Liter);
    
end

%begin loop
looplogic = true;
loopcounter = 1;
while looplogic
    
    % E-step
    
    %cycle through measurements
    for ii = 1:N
        
        %cycle through gaussians
        for jj = 1:k
            gamma(ii,jj) = w_hat(jj)*gaussEval(y(ii),mu_hat(jj),var_hat(jj));
        end
        
        %normalize
        gamma(ii,:) = gamma(ii,:)/sum(gamma(ii,:));
    end
    
    %calculate n
    n = sum(gamma,1);
    
    % M-step
    w_hat = n/N;
    for ii = 1:k
        mu_hat(ii) = sum(gamma(:,ii).*y)/n(ii);
        var_hat(ii) = sum(gamma(:,ii).*(y - mu_hat(ii)).^2)/n(ii);
    end
    
    % compute this iterations log liklihood
    L2 = 0;
    for ii = 1:N
        Liter = 0;
        for jj = 1:k
            Liter = Liter + w_hat(jj)*gaussEval(y(ii),mu_hat(jj),var_hat(jj));
        end
        L2 = L2 + log(Liter);
    end
    
    %loop management
    loopcounter = loopcounter + 1;
    if(loopcounter >= loopmax)
        looplogic = false;
        disp("Warning: EM exit from max iterations")
    end
    
    if(abs(L2 - L) <= tol)
        looplogic = false;
    else
        L = L2; %update log liklihood for next iteration
    end
    
end

%assign output
EMmodel.mu = mu_hat;
EMmodel.var = var_hat;
EMmodel.w = w_hat';
EMmodel.k = k;

end

