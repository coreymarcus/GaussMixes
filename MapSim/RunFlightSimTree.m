clear
close all
clc


%% Options

%path for sandbox
addpath("../sandbox/")
addpath("../../matlabScripts")

%how do we evaluate the GM at the end?
GMevalmeth = 'gauss'; %uses the standard gaussian
% GMevalmeth = 'mixed'; %uses a mixed representation based on s and t

%time parameters
t_sim = 59;
t_step = 1;

%trajectory
traj = 'Ramp'; %straight ramp descent

% true terrain
terrain = 'GentleParabola';

% measurement model
meas_model = 'Simple';

%initialization settings
N_meas_init = 100;
R_init = 0.1*eye(2);

%cycle settings
N_meas_cycle = 50;
FOV = pi/6;
R_cycle = 0.01*eye(2);


%% Main

%figure for plotting
% main_fig = figure;

%draw initial measurements
x1 = -400;
x2 = 400;
x_init = (x2 - x1)*rand(N_meas_init,1) + x1;
y_init = TruthEval(x_init,terrain);

%corrupt with noise
x_meas = mvnrnd(x_init, R_init(1,1)*eye(N_meas_init))';
y_meas = mvnrnd(y_init, R_init(2,2)*eye(N_meas_init))';

% Initialize a tree
tree = SlopedCell();
tree.x1_ = x1 - 10;
tree.x2_ = x2 + 10;
tree.depth_ = 0;
tree.maxDepth_ = 6;
tree.nPointsMax_ = 25;
tree = tree.AddPoints(x_init, R_init(1,1)*eye(N_meas_init), y_meas,  R_init(2,2)*eye(N_meas_init));
minwidth = (tree.x2_ - tree.x1_)*(0.5^(tree.maxDepth_));

%initialize the DEM estimate
Ndem = (tree.x2_ - tree.x1_)/minwidth;
xhatDEM = zeros(Ndem,1);
PhatDEM = ones(Ndem,1);
grid = linspace(x1-10,x2+10,Ndem);
[xhatDEM, PhatDEM] = updateDEM(xhatDEM, PhatDEM, x_meas, y_meas, grid, R_init(2,2));

%sample truth
truth_samp = [x1:.1:x2;
    TruthEval(x1:.1:x2,terrain)];

%time vector
t = 0:t_step:t_sim;
N_t = length(t);

%initialize storage
x = zeros(3,N_t);

%initialize measurement storage vector
meas = zeros(2,N_t*N_meas_cycle);

%theta values for measurement FOV
theta = linspace(-FOV/2, FOV/2, N_meas_cycle);

%Plot initial results
% tree.PlotTree(main_fig,true);

%loop
for ii = 1:N_t
    
    %get true state
    x_iter = GetVehiclePose(t(ii),traj);
    
    %measurement noise
    noise_iter = mvnrnd([0; 0],R_cycle,N_meas_cycle)';
    
    %get measurment
    meas_iter = zeros(2,N_meas_cycle);
    for jj = 1:N_meas_cycle
        meas_iter(:,jj) = GetTerrainMeasurement(x_iter,theta(jj),zeros(2),meas_model,terrain) ...
            + noise_iter(:,jj);
    end
    
    % Add measurements to the tree
    tree = tree.AddPoints(meas_iter(1,:),...
        R_cycle(1,1)*eye(N_meas_cycle),...
        meas_iter(2,:),...
        R_cycle(2,2)*eye(N_meas_cycle));
    
    % Update DEM
    [xhatDEM, PhatDEM] = updateDEM(xhatDEM, PhatDEM, meas_iter(1,:), meas_iter(2,:), grid, R_cycle(2,2));
    
    % plot
%     clf(main_fig);
%     tree.PlotTree(main_fig,true);
    
    %store
    x(:,ii) = x_iter;
    meas(:,(ii-1)*N_meas_cycle+1:ii*N_meas_cycle) = meas_iter;
    
end

%% Plotting

DEM_plot = figure;
[~, xmaxDEM, yargmaxDEM, maxDEM] = PlotDEM(DEM_plot, xhatDEM, PhatDEM, grid, x1,x2,min(truth_samp(2,:)) - 1,max(truth_samp(2,:)) + 1,.1,.1);
plot3(truth_samp(1,:),truth_samp(2,:),maxDEM*ones(length(x1:.1:x2),1),'k','LineWidth',2)
plot3(xmaxDEM, yargmaxDEM, maxDEM*ones(length(xmaxDEM),1),'r','LineWidth',2)
title('PDF from DEM')

final_fig = figure;
plot(truth_samp(1,:),truth_samp(2,:),'k','LineWidth',2)
hold on
tree.PlotTree(final_fig,false);
plot(xmaxDEM, yargmaxDEM, '--r', 'LineWidth',2)
% plot(x(1,:),x(2,:),'--r','LineWidth',1.5)

final_fig2 = figure;
plot(truth_samp(1,:),truth_samp(2,:),'k','LineWidth',2)
hold on
tree.PlotTree(final_fig2,true);
plot(xmaxDEM, yargmaxDEM, '--r', 'LineWidth',2)
% plot(x(1,:),x(2,:),'--r','LineWidth',1.5)


