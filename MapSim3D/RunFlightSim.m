clear all
close all
clc

%% Options

% seed
rng(3);

%time parameters
t0 = 0;
tf = 120;
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
somemeastarg = 5;

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

% settings to create a movie for the presentation
createmovie = true;
if(createmovie)
    movieplot = figure;
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

% Sample truth
xtruthsamp = -tree.halfWidth_:truthgsd:tree.halfWidth_;
ytruthsamp = -tree.halfHeight_:truthgsd:tree.halfHeight_;
ztruthsamp = TruthEval(xtruthsamp,ytruthsamp,terrain);

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
    %     tree.PlotTree(treeplot, 0)
    %     axis([-500 500 -500 500 -500 500])
    %     scatter3(meas_all(1,:),meas_all(2,:),meas_all(3,:),'x')
    %     pause(0)
    %     clf(treeplot)

    if(createmovie)
        % figure for tree
        figure(movieplot);
        clf(movieplot);
        tree.PlotTree(movieplot, 0)
        mesh(xtruthsamp,ytruthsamp,ztruthsamp,'FaceColor','none')
        xlabel('x [m]')
        ylabel('y [m]')
        zlabel('z [m]')
        test = axis;
        test([1 3]) = test([1 3]) - 3;
        test([2 4]) = test([2 4]) + 3;
        test(5:6) = [-15 15];
        axis(test);
        view(47,39)

        % Draw bounding box for LIDAR
        UL = meas_iter(:,1);
        LR = meas_iter(:,end);
        xbox = [UL(1), LR(1), LR(1), UL(1), UL(1)];
        ybox = [UL(2), UL(2), LR(2), LR(2), UL(2)];
%         zbox = 0.5 + [UL(3), LR(3), LR(3), UL(3), UL(3)];
        zbox = ones(1,5);
        plot3(xbox,ybox,zbox,'r','LineWidth',2)
        saveas(gcf,strcat('figs/movie/movie',num2str(ii,'%03d'),'.png'))
    end
end

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
zestplot(abs(zestplot) > 2.5) = NaN;

% Set large errors to zero
maperr2(abs(maperr2) > 20) = NaN;

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

% Find the truth surface normals
[Xtruthsamp, Ytruthsamp] = meshgrid(xtruthsamp,ytruthsamp);
[truthnormX, truthnormY, truthnormZ] = surfnorm(Xtruthsamp,...
    Ytruthsamp,ztruthsamp);

% Query estimated surface normals
surfnormgsd = truthgsd;
estsurfnorm = NaN*zeros(length(ytruthsamp),length(xtruthsamp),3);
estsurfnorm = tree.QueryMap(surfnormgsd, estsurfnorm, 'SurfNorm');

% Find error in surface normals
surfnormerr = zeros(length(ytruthsamp),length(xtruthsamp));
for ii = 1:length(xtruthsamp)
    for jj = 1:length(ytruthsamp)
        surfnormerr(jj,ii) = acosd(truthnormX(jj,ii)*estsurfnorm(jj,ii,1) + ...
            truthnormY(jj,ii)*estsurfnorm(jj,ii,2) + ...
            truthnormZ(jj,ii)*estsurfnorm(jj,ii,3));
    end
end

% Find estimate error for the tree
tree = tree.FindTruthPlane(truthgsd, terrain);

% Query Estimate Error Norm
esterrnormgsd = 1.0;
esterrnorm = NaN*zeros(round(2*tree.halfHeight_/esterrnormgsd),...
    round(2*tree.halfWidth_/esterrnormgsd));
esterrnorm = tree.QueryMap(esterrnormgsd, esterrnorm, 'EstErrNorm');

% Query Estimate Error
esterrgsd = 1.0;
esterr = NaN*zeros(round(2*tree.halfHeight_/esterrgsd),...
    round(2*tree.halfWidth_/esterrgsd),3);
esterr = tree.QueryMap(esterrgsd, esterr, 'EstErr');
esterrunwrap = [reshape(esterr(:,:,1),1,[]);
    reshape(esterr(:,:,2),1,[]);
    reshape(esterr(:,:,3),1,[])];

% Query Estimate Variance
estvargsd = 1.0;
estvar = NaN*zeros(round(2*tree.halfHeight_/estvargsd),...
    round(2*tree.halfWidth_/estvargsd),3);
estvar = tree.QueryMap(estvargsd, estvar, 'EstVar');
estvarunwrap = [reshape(estvar(:,:,1),1,[]);
    reshape(estvar(:,:,2),1,[]);
    reshape(estvar(:,:,3),1,[])];

%% Save results
filename = ['results/', datestr(datetime('now'),'yyyy-mm-dd-HH-MM')];
save(filename,'-regexp', '^(?!(movieplot)$).')

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
%plot3(trajsamp(1,:),trajsamp(2,:),trajsamp(3,:),'k-*','LineWidth',2)
scatter3(meas_all(1,:),meas_all(2,:),meas_all(3,:),'x')
legend('Truth Terrain','Spacecraft Trajectory','Location','northeast')
axis equal
saveas(gcf,'figs/truth.pdf')

% Plot truth only
truthplot2 = figure;
mesh(xtruthsamp,ytruthsamp,ztruthsamp,'FaceColor','interp')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
view(2)
cbar5 = colorbar;
cbar5.Label.Interpreter = 'latex';
cbar5.Label.String = 'Elevation [m]';
saveas(gcf,'figs/truthelev.pdf')
title('Truth')


% Plot truth with cell divisions
truthplot3 = figure;
mesh(xtruthsamp,ytruthsamp,ztruthsamp,'FaceColor','interp')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
view(2)
cbar5 = colorbar;
cbar5.Label.Interpreter = 'latex';
cbar5.Label.String = 'Elevation [m]';
hold on
tree.PlotTree(truthplot3, 100)
saveas(gcf,'figs/truthpluscells.pdf')
title('Truth With Cell Divisions')

% Plot estimate
zestplot = zestsamp;
zestplot(zestplot > 2.5) = 2.5;
zestplot(zestplot < -0.4) = -0.4;
estplot = figure;
mesh(xtruthsamp,ytruthsamp,zestplot,'FaceColor','flat')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
view(2)
cbar0 = colorbar;
cbar0.Label.Interpreter = 'latex';
cbar0.Label.String = 'Elevation [m]';
saveas(gcf,'figs/estelev.pdf')
title('Estimate')

% Plot tree and truth in x-y
% xyplot = figure;
% tree.PlotTree(xyplot, 0);
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
tree.PlotTree(treeplot, 0)
% scatter3(meas_all(1,:),meas_all(2,:),meas_all(3,:),'x')
% axis equal
mesh(xtruthsamp,ytruthsamp,ztruthsamp,'FaceColor','none')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
test = axis;
test(5:6) = [-15 15];
axis(test);
view(47,39)
saveas(gcf,'figs/Mfinal3d.pdf')


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

% Determine maperror 2 colors
% errinterp = [-7 -1 -.5:.1:.5 1 30];
errinterp = [0:.1:1 2 10 30];
colorinterp = linspace(0,1,length(errinterp));
maperr2color = zeros(size(maperr2));
for ii = 1:length(xtruthsamp)
    for jj = 1:length(ytruthsamp)
        maperr2color(jj,ii) = interp1(errinterp,colorinterp,abs(maperr2(jj,ii)));
    end
end

% Create tick labels
ticklabels = num2cell(errinterp);
for ii = 1:length(ticklabels)
    ticklabels{ii} = num2str(ticklabels{ii});
end

% Plot a more natural error
errplot2 = figure;
mesh(xtruthsamp,ytruthsamp,maperr2,maperr2color,'FaceColor','interp')
xlabel('x [m]')
ylabel('y [m]')
cbar2 = colorbar;
cbar2.Label.Interpreter = 'latex';
cbar2.Label.String = 'Absolute Error [m]';
cbar2.Ticks = colorinterp;
cbar2.TickLabels = ticklabels;
view(2)
saveas(gcf,'figs/abserr.pdf')
title('Absolute Error')

% Determine surface normal error colors
errinterp2 = [0:10 30 45 180];
colorinterp2 = linspace(0,1,length(errinterp2));
surfnormerrcolor = zeros(size(surfnormerr));
for ii = 1:length(xtruthsamp)
    for jj = 1:length(ytruthsamp)
        surfnormerrcolor(jj,ii) = interp1(errinterp2,colorinterp2,abs(surfnormerr(jj,ii)));
    end
end

% Create tick labels
ticklabels2 = num2cell(errinterp2);
for ii = 1:length(ticklabels2)
    ticklabels2{ii} = num2str(ticklabels2{ii});
end

% Plot surface normal error
surfnormplot = figure;
mesh(xtruthsamp,ytruthsamp,surfnormerr,surfnormerrcolor,'FaceColor','interp')
xlabel('x [m]')
ylabel('y [m]')
cbar7 = colorbar;
cbar7.Label.Interpreter = 'latex';
cbar7.Label.String = 'Error [deg.]';
cbar7.Ticks = colorinterp2;
cbar7.TickLabels = ticklabels2;
view(2)
saveas(gcf,'figs/surfnormerr.pdf')
title('Surface Normal Error')

% Plot fitscore
fitplot = figure;
mesh(fitscore,'FaceColor','flat')
view(2)
cbar3 = colorbar;
title('Fit Score')

% Plot estimate error norm
esterrnormplot = figure;
mesh(esterrnorm,'FaceColor','flat')
view(2)
cbar6 = colorbar;
title('Estimate Error Norm')

% Plot estimate error
esterrplot = figure;
for ii = 1:3
    subplot(3,1,ii)
    mesh(esterr(:,:,ii),'FaceColor','flat')
    hold on
    mesh(3*sqrt(estvar(:,:,ii)),"EdgeColor","r","FaceColor","none")
    mesh(-3*sqrt(estvar(:,:,ii)),"EdgeColor","r","FaceColor","none")
    view(2)
    cbar6 = colorbar;
    title(strcat('Estimate Error Element ',num2str(ii)))
end

% Plot estimate unwrapped error
esterrunwrapplot = figure;
for ii = 1:3
    subplot(3,1,ii)
    plot(esterrunwrap(ii,:),'k')
    hold on
    plot(3*sqrt(estvarunwrap(ii,:)),'r')
    plot(-3*sqrt(estvarunwrap(ii,:)),'r')
    title(strcat('Estimate Error Unwrapped Element ',num2str(ii)))
    legend("Error","3 \sigma interval")
end

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