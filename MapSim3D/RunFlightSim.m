clear all
close all
clc

%% Options

% seed
rng(3);

%time parameters
t0 = 0;
tf = 100;
t_step = 1;
pause_length = 0; % time to pause and look at plots for

%trajectory
% traj = 'Ramp'; %straight ramp descent
% traj = 'ShortRamp'; % Shorter
traj = 'Scan'; % Duplicating the ALHAT trajectory

% true terrain
% terrain = 'GentleParabola';
% terrain = 'FlatBottomParabola';
terrain = 'RocksAndCraters';

% measurement model
meas_model = 'Simple';

%sensor settings
% Nmeaselev = 50;
% Nmeasaz = 50;
% elevFOV = pi/6;
% azFOV = pi/6;
% Rmeas = .01*eye(3);

% approximate alhat sensor settings
Nmeaselev = 128;
Nmeasaz = 128;
elevFOV = pi/180;
azFOV = pi/180;
Rmeas = (0.08^2)*eye(3);

%other settings
truthgsd = .5;
somemeasonly = false; % run entire trajectory or just first measurment
somemeastarg = 1;

% tree settings
switch terrain
    case 'RocksAndCraters'
        treewidth = 60;
        treeheight = 60;
        maxdepth = 6;
    otherwise
        treewidth = 1400;
        treeheight = 1400;
        maxdepth = 4;
end


%% Main

% Time parameters
tsamp = t0:t_step:tf;
Nsamp = length(tsamp);

% Sensor parameters
elevangles = linspace(elevFOV/2,-elevFOV/2,Nmeaselev);
azangles = linspace(azFOV/2,-azFOV/2,Nmeasaz);
Nmeas = Nmeaselev*Nmeasaz;

% Get bearings in the body frame
bmatbody = zeros(3,Nmeas);
for ii = 1:Nmeasaz
    for jj = 1:Nmeaselev
        bxz = cos(elevangles(jj));
        biter = [-bxz*sin(azangles(ii)), -sin(elevangles(jj)), bxz*cos(azangles(ii))]';
        bmatbody(:,ii + (jj - 1)*Nmeasaz) = biter;
    end
end

% Initialize the tree
tree = PlanarCell();
tree.center_ = [0 0]';
tree.halfWidth_ = treewidth/2;
tree.halfHeight_ = treeheight/2;
tree.nPointsMax_ = 100;
tree.maxDepth_ = maxdepth;
tree.depth_ = 0;

% Get truth trajectory at all timesteps
trajsamp = GetVehiclePose(tsamp,traj);

% store measuremnets
meas_all = [];

% Begin the main loop
if somemeasonly
    iters = somemeastarg;
else
    iters = Nsamp - 1;
end

for ii = 1:iters
    
    % Rotate bmatbody into the inertial frame
    quatmap2body = trajsamp(4:7,ii)';
    bhatinertial = quatrotate(quatconj(quatmap2body),bmatbody')';
    
    % Get truth terrain measurements
    meas_iter = zeros(3,Nmeas);
    Rmeas_iter = zeros(3,3,Nmeas);
    m0mat = zeros(3,Nmeas);
    
    tic
    for jj = 1:Nmeas
        meas_iter(:,jj) = GetTerrainMeasurement(trajsamp(1:3,ii),...
            bhatinertial(:,jj), meas_model, terrain);
        Rmeas_iter(:,:,jj) = Rmeas;
        m0mat(:,jj) = trajsamp(1:3,ii);
    end
    toc
    
    % Get measurement noise
    meas_noise = mvnrnd([0 0 0]',Rmeas,Nmeas)';
    meas_iter = meas_iter + meas_noise;
    
    
    % Store measurements
    meas_all = [meas_all, meas_iter];
    
    % Corrupt bhat based on noisy measurement
    bhatmeas = zeros(3,Nmeas);
    relpos = meas_iter - m0mat;
    for jj = 1:Nmeas
        bhatmeas(:,jj) = relpos(:,jj)/norm(relpos(:,jj));
    end
    
    % Add measurements to the tree in a randomized fashion
    Nadd = 15;
    addidx = randi(Nadd,1,Nmeas);
    tic
    for jj = 1:Nadd
        targs = addidx == jj;
        tree = tree.AddPoints(meas_iter(:,targs),...
            Rmeas_iter(:,:,targs),...
            bhatmeas(:,targs),...
            m0mat(:,targs));
    end
    toc
    % Plot
    %     tree.PlotTree(treeplot)
    %     axis([-500 500 -500 500 -500 500])
    %     scatter3(meas_all(1,:),meas_all(2,:),meas_all(3,:),'x')
    %     pause(0)
    %     clf(treeplot)
end

% Sample truth
xtruthsamp = -tree.halfWidth_:truthgsd:tree.halfWidth_;
ytruthsamp = -tree.halfHeight_:truthgsd:tree.halfHeight_;
ztruthsamp = TruthEval(xtruthsamp,ytruthsamp,terrain);

% Sample error
zestsamp = zeros(size(ztruthsamp));
for ii = 1:length(xtruthsamp)
    for jj = 1:length(ytruthsamp)
        zestsamp(jj,ii) = tree.GetZCoordinate(xtruthsamp(ii),ytruthsamp(jj));
    end
end
maperr = log10(abs(ztruthsamp - zestsamp));
maperr2 = ztruthsamp - zestsamp;

% Trim estimate
zestplot = zestsamp;
zestplot(abs(zestplot) > 15) = NaN;

% Set large errors to zero
maperr2(abs(maperr2) > 50) = NaN;

% Query map fit score
fitscoregsd = 1.0;
fitscore = NaN*zeros(round(2*tree.halfHeight_/fitscoregsd),...
    round(2*tree.halfWidth_/fitscoregsd));
fitscore = tree.QueryMap(fitscoregsd, fitscore, 'AvgFitScore');

% Query map points
pointcountgsd = 1.0;
pointscount = NaN*zeros(round(2*tree.halfHeight_/pointcountgsd),...
    round(2*tree.halfWidth_/pointcountgsd));
pointscount = tree.QueryMap(pointcountgsd, pointscount, 'NumPoints');

% Query rejected measurements
rejectgsd = 1.0;
rejectcount = NaN*zeros(round(2*tree.halfHeight_/rejectgsd),...
    round(2*tree.halfWidth_/rejectgsd));
rejectcount = tree.QueryMap(rejectgsd, rejectcount, 'NumReject');

% Save results
filename = ['results/', datestr(datetime('now'),'yyyy-mm-dd-HH-MM')];
save(filename)

%% Plotting

set(0,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% Plot truth
truthplot = figure;
mesh(xtruthsamp,ytruthsamp,ztruthsamp,'FaceColor','flat')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Truth')
hold on
plot3(trajsamp(1,:),trajsamp(2,:),trajsamp(3,:),'k-*','LineWidth',2)
scatter3(meas_all(1,:),meas_all(2,:),meas_all(3,:),'x')
legend('Truth Terrain','Spacecraft Trajectory','Location','northeast')
axis equal
saveas(gcf,'figs/truth.pdf')

% Plot truth only
truthplot2 = figure;
mesh(xtruthsamp,ytruthsamp,ztruthsamp,'FaceColor','flat')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Truth')
view(2)
cbar5 = colorbar;

% Plot estimate
estplot = figure;
mesh(xtruthsamp,ytruthsamp,zestplot,'FaceColor','flat')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
title('Estimate')
view(2)
cbar0 = colorbar;

% Plot tree and truth in x-y
% xyplot = figure;
% tree.PlotTree(xyplot);
% plot3(trajsamp(1,:),trajsamp(2,:),trajsamp(3,:),'LineWidth',2)
% xlabel('x [m]')
% ylabel('y [m]')
% startpoint = scatter3(trajsamp(1,1),trajsamp(2,1),trajsamp(3,1),'filled');
% startpoint.MarkerFaceColor = 'r';
% startpoint.MarkerEdgeColor = 'k';
% endpoint =  scatter3(trajsamp(1,end),trajsamp(2,end),trajsamp(3,end)+1,'filled');
% endpoint.MarkerFaceColor = 'g';
% endpoint.MarkerEdgeColor = 'k';
% plot3(200*cos(0:.1:2.1*pi),200*sin(0:.1:2.1*pi),50*ones(1,length(0:.1:2.1*pi)),'LineWidth',2)
% view(2)
% axis equal
% legend('Spacecraft Trajectory','Spacecraft Start','Spacecraft Stop','Boundary of Flat Terrain')
% saveas(gcf,'figs/Mfinalxy.pdf')

% figure for tree
treeplot = figure;
tree.PlotTree(treeplot)
% scatter3(meas_all(1,:),meas_all(2,:),meas_all(3,:),'x')
% axis equal
mesh(xtruthsamp,ytruthsamp,ztruthsamp,'FaceColor','none')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
% axis([-600 600 -600 600 -10 500])
view(0,25)
saveas(gcf,'figs/Mfinal3d.pdf')
test = axis;
test(5:6) = [-15 15];
axis(test);

% Plot error
errplot = figure;
mesh(xtruthsamp,ytruthsamp,maperr,'FaceColor','interp')
xlabel('x [m]')
ylabel('y [m]')
cbar = colorbar;
cbar.Label.Interpreter = 'latex';
cbar.Label.String = '$\log_{10}(|\mathcal{M} - \mathcal{S}|)$ [m]';
view(2)
saveas(gcf,'figs/Mfinalerr.pdf')
title('Error')
% axis([-600 600 -600 600 -10 10])

% Plot a more natural error
errplot2 = figure;
mesh(xtruthsamp,ytruthsamp,maperr2,'FaceColor','interp')
xlabel('x [m]')
ylabel('y [m]')
title('Error 2')
cbar2 = colorbar;
cbar2.Label.Interpreter = 'latex';
cbar2.Label.String = '$|\mathcal{M} - \mathcal{S}$ [m]';
view(2)

% Plot fitscore
fitplot = figure;
mesh(fitscore,'FaceColor','flat')
view(2)
cbar3 = colorbar;
title('Fit Score')

% Plot number of points
ptsplot = figure;
mesh(pointscount,'FaceColor','flat')
view(2)
cbar4 = colorbar;
title('Number of Measurements')

% Plot number of rejects
rejectplot = figure;
mesh(rejectcount,'FaceColor','flat')
view(2)
cbar5 = colorbar;
title('Number of Rejected Measurements')