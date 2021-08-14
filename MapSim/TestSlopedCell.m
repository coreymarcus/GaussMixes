clear
close all
clc

%% Options
addpath("../../matlabScripts")

% Seed
% rng(3)

% truth shape
% shape = 'Parabola';
shape = 'Rock';


%% Test division

% Initialize a sloped cell
tree = SlopedCell();
tree.x1_ = -5;
tree.x2_ = 5;
tree.depth_ = 0;
tree.maxDepth_ = 3;

% Assign some random points to the cell
var_x = .001;
var_y = .001;
tree.points_(1,:) = 10*(rand(1,25) - 0.5);
tree.points_(2,:) = TruthEval(tree.points_(1,:), shape) + mvnrnd(0,var_y,25)';
tree.points_(1,:) = tree.points_(1,:) + mvnrnd(0,var_x,25)';

% Divide the tree
tree = tree.Divide();

% Plot to see how we did
fighandle = figure;
tree.PlotTree(fighandle, true);

%% Test Update

% Params
Nadd = 10;
Ncycle = 15;

for ii = 1:Ncycle
    
    %Create measurements
    pts2add = zeros(2,Nadd);
    pts2add(1,:) = 10*(rand(1,Nadd) - 0.5);
    pts2add(2,:) = TruthEval(pts2add(1,:),shape) + mvnrnd(0,var_y,Nadd)';
    pts2add(1,:) = pts2add(1,:) + mvnrnd(0,var_x,Nadd)';
    
    % Add
    tree = tree.AddPoints(pts2add(1,:), var_x*eye(Nadd), pts2add(2,:), var_y*eye(Nadd));
end



% Plot
fighandle2 = figure;
tree.PlotTree(fighandle2, true);