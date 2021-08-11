function [plane_hat, cov] = PlaneEstimateML(r_mat,cov_r_mat)
%PlaneEstimateML

% locals
npts = size(r_mat,2);
% cov_eta_mat = zeros(npts,1);
% A = zeros(npts,1);
% D = zeros(3,npts);
J = [eye(2); zeros(1,2)];

% Initialize estimate
nhat = [0 0 1]';
d = 0;

% Convergence criteria
max_iters = 100;
tol = 1E-9;

% loop
for ii = 1:max_iters
    
    % find nhatplus
    kappa = nhat(1:2);
    lambda = nhat(3);
    nhatplus = zeros(3);
    nhatplus(1:2,1:2) = eye(2) - kappa*kappa'/(1+lambda);
    nhatplus(3,1:2) = -kappa';
    nhatplus(1:2,3) = kappa;
    nhatplus(3,3) = lambda;
    
    %Find some metrics for each measurement
    b = zeros(3,1);
    H = zeros(3,3);
    for jj = 1:npts
        cov_eta = nhat'*cov_r_mat(:,:,jj)*nhat;
        A = r_mat(:,jj)'*nhat - d;
        D = [-r_mat(:,jj)'*nhatplus*J, 1];
        b = b + D'*A/cov_eta;
        H = H + D'*D/cov_eta;
    end
    
    % find optimal step
    step = H\b;
    
    % update
    phi = norm(step(1:2));
    if(phi == 0)
        expphi = [step(1:2); 1];
    else
        expphi = [sin(phi)*step(1:2)/phi;
            cos(phi)];
    end
    nhat = nhatplus*expphi;
    d = d + step(3);
    
    % check for convergence
    if(norm(step) < tol)
        break;
    end
    
end

plane_hat = [nhat; d];
cov = inv(H);

end

