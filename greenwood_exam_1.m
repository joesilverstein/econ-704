clear;
close all;

%% 4 %%

alpha = 1/3;
phi = 0.5;
tau = 0.2;
Z = @(h) (1-alpha)*phi*h^(-1) - 2*(1-alpha)^2*phi*tau*h^(-alpha) - (1-phi)/(1-h);
h = bisection(Z)
w = (1-alpha)*h^(-alpha)
r = alpha*h^(1-alpha)
y = h^(1-alpha)
lambda = tau*(w*h)^2

%% 5 %%

N = 101;
tau = linspace(0,2,N);
h5 = zeros(N,1); y5 = h5; lambda5 = h5;
for i=1:N
    Z = @(h) (1-alpha)*phi*h^(-1) - 2*(1-alpha)^2*phi*tau(i)*h^(-alpha) - (1-phi)/(1-h);
    h5(i) = bisection(Z); % labor supply
    y5(i) = h5(i)^(1-alpha); % total production
    w = (1-alpha)*h5(i)^(-alpha);
    lambda5(i) = tau(i)*(w*h5(i))^2; % government revenue
end;

%% 6 %%

N = 101;
tau = linspace(0,2,N);
h6 = zeros(N,1); y6 = h6; lambda6 = h6;
for i=1:N
    Z = @(h) ((1-alpha)*phi*h^(-alpha)-2*(1-alpha)^2*phi*tau(i)*h^(1-2*alpha)) / (h^(1-alpha)-(1-alpha)^2*tau(i)*h^(2*(1-alpha))) - (1-phi) / (1-h);
    h6(i) = bisection(Z); % labor supply
    y6(i) = h6(i)^(1-alpha); % total production
    w = (1-alpha)*h6(i)^(-alpha);
    lambda6(i) = tau(i)*(w*h6(i))^2; % government revenue
end;
plot(tau, h5, tau, h6, tau, y5, tau, y6, tau, lambda5, tau, lambda6)
title('Labor Supply, Total Production, and Government Revenue as Functions of \tau')
legend('Labor Supply (Rebate)','Labor Supply (No Rebate)','Total Production (Rebate)','Total Production (No Rebate)','Government Revenue (Rebate)','Government Revenue (No Rebate)')
