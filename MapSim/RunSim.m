%Runs the mapping simulation
clear
close all
clc


%% Options

%path for sandbox
addpath("../sandbox/")


%% Main

%generate some initial points
x1 = -5;
x2 = 5;
Ninit = 50;
x_init = (x2 - x1)*rand(Ninit,1) + x1;
y_init = TruthEval(x_init);

%corrupt with noise
sig2 = 0.1;
x_meas = mvnrnd(x_init, sig2*eye(Ninit))';
y_meas = mvnrnd(y_init, sig2*eye(Ninit))';

%use k-means to initialize the map
[mu_init, P_init, w_init, group_idx] = kMeans(2, [x_meas, y_meas]');

% P = GaussElement(10);
% 
% P = P.method1(4)