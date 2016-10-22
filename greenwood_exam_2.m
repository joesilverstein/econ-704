% Greenwood Take-Home Exam 2
% Joseph Silverstein

cd('C:\Users\Joe\Google Drive\Econ 704')

rho = 2;
alpha = 0.3;
delta = 0.1;
beta = 0.96;
theta = 0.6;
tauh = 0.2;
tauk = 0.2;

rss = (1-beta)/(beta*(1-tauk))+delta;
wss = (1-alpha)*((1/alpha)*rss)^(alpha/(alpha-1));
hss = ((1-tauh)*wss)^(1/theta);
kss = ((1/alpha)*rss)^(1/(alpha-1))*hss;
yss = kss^alpha*hss^(1-alpha);

k1 = kss/2;

% We want to solve the boundary value problem starting at k1 and ending 
% at kT by iterating backwards. To do this, we only need to determine the
% optimal kT1 (i.e., k_T-1).

% I solved for h using pencil and paper. Here is the formula:
h = @(k) ((1-tauh)*(1-alpha)*k^alpha)^(1/(theta+alpha));

r = @(k) alpha*k^(alpha-1)*h(k)^(1-alpha);
w = @(k) (1-alpha)*k^alpha*h(k)^(-alpha);

lambda = @(k) tauh*w(k)*h(k) + tauk*(r(k)-delta)*k;

E = @(k,kp,kpp) ((1-tauh)*w(k)*h(k) + (1-tauk)*(r(k)-delta)*k + k + lambda(k) - kp - h(k)^(1+theta)/(1+theta))^(-rho) - beta * ((1-tauk)*(r(kp)-delta)+1) * ((1-tauh)*w(kp)*h(kp) + (1-tauk)*(r(kp)-delta)*kp + kp + lambda(kp) - kpp - h(kp)^(1+theta)/(1+theta))^(-rho); % Euler eqn

%% 4 %%

T = 200; % initial guess of the number of periods T required for the transition path
kT = kss;
kT1 = kss; % initial guess of k_T-1
Z = @(k) k1 - K1(k,kT,T,E); % Z defined to be strictly decreasing
kT1 = fzero(Z,kT1); % solve for K_T-1 starting from initial guess

kPath = path(kT1,kT,T,E);
hPath = zeros(T,1);
yPath = zeros(T,1);
cPath = zeros(T,1);
for t=1:T-1
    hPath(t) = h(kPath(t));
    yPath(t) = kPath(t)^alpha*hPath(t)^(1-alpha);
    cPath(t) = yPath(t) + (1-delta)*kPath(t) - kPath(t+1);
end
hPath(T) = h(kPath(T));
yPath(T) = kPath(T)^alpha*hPath(T)^(1-alpha);
cPath(T) = yPath(T) + (1-delta)*kPath(T) - kss;

figure(1)
X = 1:T;
plot(X,kPath,X,hPath,X,yPath,X,cPath)
legend('k','h','y','c')
title('Transition Paths')
xlabel('Period')
ylabel('Level')

%% 5 %%

%% Congressman K's Plan:

taukK = 0.1;

% Now we start in the old steady state and transition into the new steady
% state.
k1K = kss; % initial guess of k_T-1

rssK = (1-beta)/(beta*(1-taukK))+delta;
wssK = (1-alpha)*((1/alpha)*rssK)^(alpha/(alpha-1));
hssK = ((1-tauh)*wssK)^(1/theta);
kssK = ((1/alpha)*rssK)^(1/(alpha-1))*hssK;

kTK = kssK; % ends in the new steady state

lambdaK = @(k) tauh*w(k)*h(k) + tauk*(r(k)-delta)*k;
E = @(k,kp,kpp) ((1-tauh)*w(k)*h(k) + (1-taukK)*(r(k)-delta)*k + k + lambdaK(k) - kp - h(k)^(1+theta)/(1+theta))^(-rho) - beta * ((1-taukK)*(r(kp)-delta)+1) * ((1-tauh)*w(kp)*h(kp) + (1-taukK)*(r(kp)-delta)*kp + kp + lambdaK(kp) - kpp - h(kp)^(1+theta)/(1+theta))^(-rho); % Euler eqn

ZK = @(k) k1K - K1(k,kTK,T,E);
kT1K = fzero(ZK,kssK); % solve for K_T-1 starting from initial guess of the new steady state

kPathK = path(kT1K,kTK,T,E);
hPathK = zeros(T,1);
yPathK = zeros(T,1);
cPathK = zeros(T,1);
for t=1:T-1
    hPathK(t) = h(kPathK(t));
    yPathK(t) = kPathK(t)^alpha*hPathK(t)^(1-alpha);
    cPathK(t) = yPathK(t) + (1-delta)*kPathK(t) - kPathK(t+1);
end
hPathK(T) = h(kPathK(T));
yPathK(T) = kPathK(T)^alpha*hPathK(T)^(1-alpha);
cPathK(T) = yPathK(T) + (1-delta)*kPathK(T) - kssK;

%% Congressman L's Plan:

tauhL = 0.175;

% Now we start in the old steady state and transition into the new steady
% state.
k1L = kss; % initial guess of k_T-1

wssL = (1-alpha)*((1/alpha)*rss)^(alpha/(alpha-1));
hssL = ((1-tauhL)*wssL)^(1/theta);
kssL = ((1/alpha)*rss)^(1/(alpha-1))*hssL;

kTL = kssL; % ends in the new steady state

hL = @(k) ((1-tauhL)*(1-alpha)*k^alpha)^(1/(theta+alpha));
rL = @(k) alpha*k^(alpha-1)*hL(k)^(1-alpha);
wL = @(k) (1-alpha)*k^alpha*hL(k)^(-alpha);
lambdaL = @(k) tauh*wL(k)*hL(k) + tauk*(rL(k)-delta)*k;
E = @(k,kp,kpp) ((1-tauhL)*wL(k)*hL(k) + (1-tauk)*(rL(k)-delta)*k + k + lambdaL(k) - kp - hL(k)^(1+theta)/(1+theta))^(-rho) - beta * ((1-tauk)*(rL(kp)-delta)+1) * ((1-tauhL)*wL(kp)*hL(kp) + (1-tauk)*(rL(kp)-delta)*kp + kp + lambdaL(kp) - kpp - hL(kp)^(1+theta)/(1+theta))^(-rho); % Euler eqn

ZL = @(k) k1L - K1(k,kTL,T,E);
kT1L = fzero(ZL,kssL); % solve for K_T-1 starting from initial guess of the new steady state

kPathL = path(kT1L,kTL,T,E);
hPathL = zeros(T,1);
yPathL = zeros(T,1);
cPathL = zeros(T,1);
for t=1:T-1
    hPathL(t) = hL(kPathL(t));
    yPathL(t) = kPathL(t)^alpha*hPathL(t)^(1-alpha);
    cPathL(t) = yPathL(t) + (1-delta)*kPathL(t) - kPathL(t+1);
end
hPathL(T) = hL(kPathL(T));
yPathL(T) = kPathL(T)^alpha*hPathL(T)^(1-alpha);
cPathL(T) = yPathL(T) + (1-delta)*kPathL(T) - kssL;

%% Plot both transition paths on the same graph:

figure(2)
plot(X,yPathK,X,yPathL)
legend('Congressman K','Congressman L')
title('Transition Paths of GDP')
xlabel('Period')
ylabel('GDP')

figure(3)
plot(X,kPathK,X,kPathL)
legend('Congressman K','Congressman L')
title('Transition Paths of Capital')
xlabel('Period')
ylabel('Capital')

%% Calculate lifetime utility under the two plans:

% period utility
u = @(c,h) beta * (1/(1-rho)) * (c - h^(1+theta)/(1+theta))^(1-rho);

% Utility of the representative agent under K's plan:
sum = 0;
for t=1:T
    sum = sum + beta^t * u(cPathK(t),hPathK(t));
end
ussK = u(cPathK(T),hPathK(T)); % steady state period utility
% add utility up to 10000 periods
for t=T+1:10000
    sum = sum + beta^t * ussK;
end
UK = sum; % lifetime utility

% Utility of the representative agent under L's plan:
sum = 0;
for t=1:T
    sum = sum + beta^t * u(cPathL(t),hPathL(t));
end
ussL = u(cPathL(T),hPathL(T)); % steady state period utility
% add utility up to 10000 periods
for t=T+1:10000
    sum = sum + beta^t * ussL;
end
UL = sum;

