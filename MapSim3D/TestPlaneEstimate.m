clear
close all
clc

%% Settings

nMC = 1000;

nMeas = 100;
Rmeas = eye(3);
plane_truth = [0 0 1 1]';
plane_truth(1:3) = plane_truth(1:3)/norm(plane_truth(1:3));

xbounds = [-5 5];
ybounds = [-5 5];


%% Main

% MC parameters
cov_mat = zeros(3,3,nMC);
est_mat = zeros(4,nMC);
e_mat = zeros(3,nMC);

% MC loop
for jj = 1:nMC
    %draw sample points
    pts = zeros(3,nMeas);
    ptsCov = zeros(3,3,nMeas);
    xtrue = (xbounds(2) - xbounds(1))*rand(nMeas,1) + xbounds(1);
    ytrue = (ybounds(2) - ybounds(1))*rand(nMeas,1) + ybounds(1);
    ztrue = -(plane_truth(1)*xtrue + plane_truth(2)*ytrue - plane_truth(4))/plane_truth(3);
    rtrue = [xtrue, ytrue, ztrue]';
    r = zeros(3,nMeas);
    rCov = zeros(3,3,nMeas);
    for ii = 1:nMeas
        r(:,ii) = mvnrnd(rtrue(:,ii),Rmeas);
        rCov(:,:,ii) = Rmeas;
    end
    
    %call algorithm
    [plane_hat, cov_hat] = PlaneEstimateML(r,rCov);
    
    %store
    cov_mat(:,:,jj) = cov_hat;
    est_mat(:,jj) = plane_hat;
    e_mat(1:2,jj) = UnitVectorSubtract(plane_hat(1:3), plane_truth(1:3));
    e_mat(3,jj) = plane_hat(4) - plane_truth(4);
    
end

% See how well estimator does
eMean = mean(e_mat,2)
eCov = cov(e_mat')
predCov = mean(cov_mat,3)

%% Plot

figure
scatter3(xtrue,ytrue,ztrue)
axis equal