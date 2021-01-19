function [MAP_model] = EM_MAP(prior_model, y, options)
%EM_MAP performs a maximum a posteriori EM update given some measurements
%and a prior on the model. Sourced from "Bayesian regularization for normal
%mixture estimation and model-based clustering" [2007]

%extract local variables
k = prior_model.k; %number of gaussians
mu_hat = prior_model.mu; %means
var_hat = prior_model.var; %variances
w_hat = prior_model.w; %weights
N = length(y); %number of measurements
gamma = zeros(N,k); %prob that (row) measurement from (col) gaussian

%initialize prior parameters
mu_p = mean(y);
kappa_p = 0.01;
nu_p = 3;
zeta_p2 = var(y)/(k^2);

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
        y_bar = sum(gamma(:,ii).*y)/n(ii);
        mu_hat(ii) = (n(ii)*y_bar + kappa_p*mu_p)/(kappa_p + n(ii));
        term2 = kappa_p*n(ii)*(y_bar - mu_p)^2/(kappa_p+n(ii));
        term3 = sum(gamma(:,ii).*(y - y_bar).^2);
        var_hat(ii) = (zeta_p2 + term2 + term3)/(nu_p + n(ii) + 3);
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
MAP_model.mu = mu_hat;
MAP_model.var = var_hat;
MAP_model.w = w_hat';
MAP_model.k = k;


end

