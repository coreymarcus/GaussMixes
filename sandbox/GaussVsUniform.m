clear
close all
clc

%path
addpath('../matlabScripts/');

%bounds
x1 = 1;
x2 = 100;

%mean and var
var = (x2 - x1)^2/12;
mu = (x2 + x1)/2;

%sample
Nsample = 50;
x = linspace(x1,x2,Nsample);
unifdist = zeros(Nsample,1);
normdist = zeros(Nsample,1);
for ii = 1:Nsample
    unifdist(ii) = 1/(x2 - x1);
    normdist(ii) = gaussEval(x(ii), mu, var);
end

%plot
figure
plot(x,unifdist,x,normdist)
axis([x1 x2 0 1.5*max(normdist)])