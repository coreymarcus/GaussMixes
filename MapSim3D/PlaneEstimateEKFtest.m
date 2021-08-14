clear
close all
clc

% Create the truth plane
nhat_true = [.1 .1 1]';
nhat_true = nhat_true/norm(nhat_true);
d_true = 2;

%seed
rng(1)

% Observer location
m0 = [0 0 10]';
P_m0 = 10*eye(3); % covariance for peturbing m0

% Bearing parameters
bhat = [0 0 -1]';
P_bhat = .2*eye(3); % covariance for peturbing b

% measurement parameters
P_r = .1*eye(3); % measurement covariance

% Create a bunch of initial measurements
Ninit = 10;
rinit = zeros(3,Ninit);
Pinit = zeros(3,3,Ninit);
bpeturb = mvnrnd([0 0 0]',P_bhat,Ninit)';
m0peturb = mvnrnd([0 0 0]',P_m0,Ninit)';
rpeturb = mvnrnd([0 0 0]',P_r,Ninit)';
for ii = 1:Ninit
    
    % Find bhat
    bhat_iter = bhat + bpeturb(:,ii);
    bhat_iter = bhat_iter/norm(bhat_iter);
    
    m0_iter = m0 + m0peturb(:,ii);
    
    % Find the true measurement
    t_iter = (d_true - nhat_true'*m0_iter)/(nhat_true'*bhat_iter);
    rinit(:,ii) = m0_iter + t_iter*bhat_iter + rpeturb(:,ii);
    
    % Covariance
    Pinit(:,:,ii) = P_r;
    
end

% Initialize plane estimate
[xhat0, Phat0] = PlaneEstimateML(rinit,Pinit);

% Initialize recursive estimate
Nupdate = 5000;
xhat = [xhat0 zeros(4,Nupdate)];
Phat = zeros(3,3,Nupdate);
Phat(:,:,1) = Phat0;

% Begin the loop
bpeturb = mvnrnd([0 0 0]',P_bhat,Nupdate)';
rpeturb = mvnrnd([0 0 0]',P_r,Nupdate)';
m0peturb = mvnrnd([0 0 0]',P_m0,Nupdate)';
nSkipped = 0;
for ii = 1:Nupdate
    
    % Find bhat
    bhat_iter = bhat + bpeturb(:,ii);
    bhat_iter = bhat_iter/norm(bhat_iter);
    
    % peturb m0
    m0_iter = m0 + m0peturb(:,ii);
    
    % Find the true measurement
    t_iter = (d_true - nhat_true'*m0_iter)/(nhat_true'*bhat_iter);
    r_iter = m0_iter + t_iter*bhat_iter + rpeturb(:,ii);
    
    % Find the expected measurement
    nhat = xhat(1:3,ii);
    dhat = xhat(4,ii);
    t_exp = (dhat - nhat'*m0_iter)/(nhat'*bhat_iter);
    r_exp = m0_iter + t_exp*bhat_iter;
    
    % Skip if t_exp too big
    theta = acos(nhat'*bhat_iter)*180/pi;
    %if(abs(t_exp) > 100)
    if(theta < 110 || t_exp <= 0)
        xhat(1:3,ii+1) = nhat;
        xhat(4,ii+1) = dhat;
        Phat(:,:,ii+1) = Phat(:,:,ii);
        nSkipped = nSkipped + 1;
        disp(theta)
        disp(t_exp)
        disp(' ')
        continue;
    end
    
    % Find the measurement jacobian
    H = JacobianEval(bhat_iter(1),bhat_iter(2),bhat_iter(3),...
        dhat,...
        m0_iter(1),m0_iter(2),m0_iter(3),...
        nhat(1),nhat(2),nhat(3));
    
    % Find the kalman gain
    Pbar = Phat(:,:,ii);
    K = Pbar*H'/(H*Pbar*H' + P_r);
    
    % Find the update
    update = K*(r_iter - r_exp);
    
    % Perform the update
    xhat(1:3,ii+1) = UnitVectorAdd(nhat,update(1:2));
    xhat(4,ii+1) = dhat + update(3);
    Phat(:,:,ii+1) = (eye(3) - K*H)*Pbar*(eye(3) - K*H)' + K*P_r*K';
    
    %     if(ii == 2960)
    %         disp('here')
    %     end
    
end


figure
subplot(2,2,1)
plot(xhat(1,:) - nhat_true(1))

subplot(2,2,2)
plot(xhat(2,:) - nhat_true(2))

subplot(2,2,3)
plot(xhat(3,:) - nhat_true(3))

subplot(2,2,4)
plot(xhat(4,:) - d_true)

Phat(:,:,end)

nSkipped