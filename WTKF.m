function [xhat_new, Phat_new] = WTKF(y, xhat, Phat, Q, H)
%Weighted Total Kalman Filter - Based on algorithm from "A general Weighted
%Total Kalman Filter algorithm with numerical evaluation"

%Step One: Initialize
m = length(y);
n = length(xhat);
mu_rhat = zeros(n,1);
M = [eye(m*n, m*n), zeros(m*n, m)];
N = [zeros(m, n*m), eye(m)];
Im = eye(m,m);
A_i = H; %simple redefinition in order to bridge aerospace lingo with lingo in paper
Phi_i = eye(n,n); %static problem
theta_i = zeros(n,n); %no process noise since problem is static
E_rhatAi = zeros(m, n); %I believe this is zero because our predicted peturbation is zero

%Step Two: Iteratively Compute lambda_tilde, E_rhat, mu_rhat
loopcontinue = true;
maxloop = 100;
loopidx = 1;
tol = 1E-16;
lamb_tilde = zeros(m,1);
while loopcontinue
    
    %find lambda
    term1 = N - kron((mu_rhat + xhat)',Im)*M;
    lamb_iter = inv(term1*Q*term1' + A_i*(theta_i + Phi_i*Phat*Phi_i')*(A_i - E_rhatAi)')*(y - A_i*xhat);
    E_rhat = Q*term1'*lamb_iter;
    mu_rhat = (theta_i + Phi_i*Phat*Phi_i')*(A_i - E_rhatAi)'*lamb_iter;
    
    %update E_rhatAi
    vecE_rhatAi = M*E_rhat;
    E_rhatAi = reshape(vecE_rhatAi,m,n);
    
    
    %loop managment
    loopidx = loopidx + 1;
    if(loopidx > maxloop)
        disp("Warning: max loop achieved")
        loopcontinue = false;
    end
    
    if(norm(lamb_iter - lamb_tilde) <= tol)
        loopcontinue = false;
    end
    
    %update lambda
    lamb_tilde = lamb_iter;
end


%Step Three: Calculate the update
xhat_new = xhat + mu_rhat;
Phat_new = kron(theta_i + Phi_i*Phat*Phi_i',lamb_tilde')*M*Q*M'*...
    kron(theta_i + Phi_i*Phat*Phi_i',lamb_tilde')' + Phi_i*Phat*Phi_i';

end

