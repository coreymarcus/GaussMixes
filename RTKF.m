function [xhat_out, Phat_out] = RTKF(xhat, Phat, Phi, Sigma_w, Sigma_e, A, L, C_A, C_L)
%RTKF - "Robust Total Kalman Filter" Sourced from "A robust total Kalman
%filter algorithm with numerical evaluation"

%local variables
n = length(xhat);
m = length(L);
outertol = 1E-4;
innertol = 1E-4;

%propagate estimate
xbar = Phi*xhat;
Pbar = Phi*Phat*Phi' + Sigma_w;

%loop management
outeridx = 1;
outerloopmax = 100;
outerlogic = true;

%begin loop
while outerlogic
    
    %inner loop logic
    inneridx = 1;
    innerloopmax = 100;
    innerlogic = true;
    
    %initialize loop variables
    xhat_last = xbar;
    
    %begin inner loop
    while innerlogic
        G = C_L - kron(xhat_last',eye(n))*C_A;
        M = G*Sigma_e*G';
        e = Sigma_e*G'*inv(M)*(L - A*xhat_last);
        E_A = [e zeros(2,1)];
        A_tilde = A - E_A;
        L_tilde = L - E_A*xhat_last;
        K = inv(Pbar)*A_tilde'*inv(A*Pbar*A' + M);
        xhat_new = xbar + K*(L_tilde - A_tilde*xbar);
        Phat_new = (eye(n) - K*A_tilde)*Pbar*(eye(n) - K*A_tilde)' + K*M*K';
        
        %enforce Phat > 0
        Phat_new = 0.5*(Phat_new + Phat_new');
        
        %change
        delta2 = norm(xhat_new - xhat_last);
        xhat_last = xhat_new;
        
        %loop management
        inneridx = inneridx + 1;
        if(inneridx > innerloopmax)
            disp("Warning: Max inner loop!")
            innerlogic = false;
        elseif(delta2 <= innertol)
            innerlogic = false;
        end
    end
    
    %update
    Sigma_e = Sigma_e; %FIX THIS
    delta1 = 0;
    
    %loop managment
    outeridx = outeridx + 1;
    if(outeridx > outerloopmax)
        disp("Warning: Max outer loop!")
        outerlogic = false;
    elseif(delta1 <= outertol)
        outerlogic = false;
    end
    
    
end

end

