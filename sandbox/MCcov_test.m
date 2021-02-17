clear
close all
clc


%parameters
mu_x = 4;
var_x = 1;
mu_z = 2;
var_z = 0.5;
N = 1000000;
mu_w = 0;
var_w = 0.1;

%draw
x = mvnrnd(mu_x,var_x,N);
z = mvnrnd(mu_z,var_z,N);
w = mvnrnd(mu_w,var_w,N);

%measurement
y = x + x.*z + w;

%calculate mean and variance
mean(x.*(x.*z))
mean(x.*(x.*z) - mean(x.*(x.*z)))

%Pyy
var(y)
var_x + var_w + var_x*var_z + var_z*mu_x^2

%test
q = x.*z;

figure
hist(q,100)
