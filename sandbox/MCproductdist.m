% Lets do a little monte carlo to test our product distribution

clear
close all
clc

%Setup
N = 100000; %number of MC
mu1 = 1; %average of non-zero RV
mu2 = 2; %average of non-zero RV
var1 = .1; %variance
var2 = .2; %variance 2
n_max = 50; %this controls the number of terms we consider in our pdf

%combinatorial binomial coeff anonymous function
bicoeff = @(n,k) factorial(n)/(factorial(k)*factorial(n-k));

%pdf eval points
% x_eval = (-3+(mu1*mu2)):.01:(4+(mu1*mu2));
x_eval = -2:.01:6;

%remove x = 0
x_eval = x_eval(find(x_eval));

L = length(x_eval);
pdf_eval = zeros(1,L);
parfor ii = 1:L
    
    %initialize
    sum_eval = 0;
    
    %loops
    for n = 0:n_max
        for m = 0:(2*n)
            
            %terms of ugly fraction
            numer = mu1^m * mu2^(2*n-m) * x_eval(ii)^(2*n-m) * abs(x_eval(ii))^(m-n);
            denomer = pi * factorial(2*n) * sqrt(var1)^(n+m+1) * sqrt(var2)^(3*n - m + 1);
            
            %evaluate sum
            sum_eval = sum_eval + bicoeff(2*n,m)*(numer/denomer)*besselk(m-n,abs(x_eval(ii))/sqrt(var1*var2));          
        end
    end
    
    %assign
    pdf_eval(ii) = exp(-0.5*(mu1^2/var1 + mu2^2/var2))*sum_eval;
end

%approximate the mean and variance
dx = x_eval(2)-x_eval(1);
mu_y = trapz(x_eval,x_eval.*pdf_eval)
var_y = trapz(x_eval,((x_eval - mu1*mu2).^2).*pdf_eval)

%sample
x1 = mvnrnd(mu1,var1,N);
x2 = mvnrnd(mu2,var2,N);
y = x1.*x2;
mean(y)
var(y)



%plotting
figure
histogram(y,N/100,'Normalization','pdf');
hold on
plot(x_eval,pdf_eval,'r','LineWidth',2)
xlabel('y = x1*x2')
ylabel('p(y)')
legend('MC','Analytic')
