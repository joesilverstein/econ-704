% Greenwood Take-Home Exam 3
% Joseph Silverstein

cd('C:\Users\Joe\Google Drive\Econ 704')

% parameters
rho = 2;
beta = 0.96;
alpha = 0.3;
delta = 0.1;
dz = 0.03;
z1 = 1 + dz;
z2 = 1 - dz;
p11 = 0.7;
p12 = 1 - p11;
p21 = p12;
p22 = 1 - p21;
TOL = 10e-4;
N = 20; % number of grid points
T = 10000;
burn = 3500; % number of points to burn in the Markov chain
gridk = linspace(2.5,3.3,N); % grid of capital

%% 1 %%
kss = ((1/alpha)*((1/beta)-1+delta))^(1/(alpha-1));

%% 3(a-c), 4 %%

[vfun1,vfun2,K1,K2,kSim,busStats] = everything(dz,p11,N,TOL,alpha,beta,delta,rho,kss,burn,T,gridk);

% plot the value function
figure(1)
plot(gridk,vfun1,gridk,vfun2)
legend('V(k,z_1)','V(k,z_2)')
title('Value Function')

% plot the policy function
figure(2)
plot(gridk,K1,gridk,K2,gridk,gridk)
legend('K(k,z_1)','K(k,z_2)')
title('Policy Function')

% distribution of capital in the steady state
figure(3)
hist(kSim(burn+1:T),50)
title('Distribution of Capital in the Stochastic Steady State')

% business cycle statistics
disp(busStats)

%% 4 %%

nz = 8;
npi = 8;

griddz = linspace(0.01,0.2,nz);
gridp11 = linspace(0.1,0.9,npi);

yStd = zeros(nz,npi); cStd = zeros(nz,npi); iStd = zeros(nz,npi);
cCorr = zeros(nz,npi); iCorr = zeros(nz,npi);
yAcorr = zeros(nz,npi); cAcorr = zeros(nz,npi); iAcorr = zeros(nz,npi);

N = 5;
T = 3000;
burn = 2000;
tic
for i=1:nz
    for j=1:npi
        [vfun1,vfun2,K1,K2,kSim,busStats] = everything(griddz(i),gridp11(j),N,TOL,alpha,beta,delta,rho,kss,burn,T,gridk);
        yStd(i,j) = busStats(1,1); cStd(i,j) = busStats(2,1); iStd(i,j) = busStats(3,1);
        cCorr(i,j) = busStats(2,2); iCorr(i,j) = busStats(3,2);
        yAcorr(i,j) = busStats(1,3); cAcorr(i,j) = busStats(2,3); iAcorr(i,j) = busStats(3,3);
    end
end
toc

% determine which business cycle statistics table looks most like the U.S.
% data:
min = 10000;
for i=1:nz
    for j=1:npi
        % minimizes the Euclidean norm of the percentage deviations from
        % the target values
        tmp = sqrt((yStd(i,j)/3.5-1)^2 + (cStd(i,j)/2.2-1)^2 + (iStd(i,j)/10.5-1)^2 + ...
            + (cCorr(i,j)/0.74-1)^2 + (iCorr(i,j)/0.68-1)^2 + ...
            + (yAcorr(i,j)/0.66-1)^2 + (cAcorr(i,j)/0.72-1)^2 + (iAcorr(i,j)/0.25-1)^2);
        if tmp < min
            iMin = i; jMin = j;
            min = tmp;
        end
    end
end

busStatsMin = zeros(3);
busStatsMin(1,1) = yStd(iMin,jMin);
busStatsMin(2,1) = cStd(iMin,jMin);
busStatsMin(3,1) = iStd(iMin,jMin);
busStatsMin(1,2) = 1;
busStatsMin(2,2) = cCorr(iMin,jMin);
busStatsMin(3,2) = iCorr(iMin,jMin);
busStatsMin(1,3) = yAcorr(iMin,jMin);
busStatsMin(2,3) = cAcorr(iMin,jMin);
busStatsMin(3,3) = iAcorr(iMin,jMin);

griddz(iMin)
gridp11(jMin)
disp(busStatsMin)

% The above gives us dz=0.037142857142857 and p11=0.671428571428571. Now we
% run again in a more local area:

nz = 8;
npi = 8;

griddz = linspace(0.025,0.045,nz);
gridp11 = linspace(0.55,0.75,npi);

yStd = zeros(nz,npi); cStd = zeros(nz,npi); iStd = zeros(nz,npi);
cCorr = zeros(nz,npi); iCorr = zeros(nz,npi);
yAcorr = zeros(nz,npi); cAcorr = zeros(nz,npi); iAcorr = zeros(nz,npi);

N = 5;
T = 3000;
burn = 2000;
tic
for i=1:nz
    for j=1:npi
        [vfun1,vfun2,K1,K2,kSim,busStats] = everything(griddz(i),gridp11(j),N,TOL,alpha,beta,delta,rho,kss,burn,T,gridk);
        yStd(i,j) = busStats(1,1); cStd(i,j) = busStats(2,1); iStd(i,j) = busStats(3,1);
        cCorr(i,j) = busStats(2,2); iCorr(i,j) = busStats(3,2);
        yAcorr(i,j) = busStats(1,3); cAcorr(i,j) = busStats(2,3); iAcorr(i,j) = busStats(3,3);
    end
end
toc

% determine which business cycle statistics table looks most like the U.S.
% data:
min = 10000;
for i=1:nz
    for j=1:npi
        % minimizes the Euclidean norm of the percentage deviations from
        % the target values
        tmp = sqrt((yStd(i,j)/3.5-1)^2 + (cStd(i,j)/2.2-1)^2 + (iStd(i,j)/10.5-1)^2 + ...
            + (cCorr(i,j)/0.74-1)^2 + (iCorr(i,j)/0.68-1)^2 + ...
            + (yAcorr(i,j)/0.66-1)^2 + (cAcorr(i,j)/0.72-1)^2 + (iAcorr(i,j)/0.25-1)^2);
        if tmp < min
            iMin = i; jMin = j;
            min = tmp;
        end
    end
end

busStatsMin = zeros(3);
busStatsMin(1,1) = yStd(iMin,jMin);
busStatsMin(2,1) = cStd(iMin,jMin);
busStatsMin(3,1) = iStd(iMin,jMin);
busStatsMin(1,2) = 1;
busStatsMin(2,2) = cCorr(iMin,jMin);
busStatsMin(3,2) = iCorr(iMin,jMin);
busStatsMin(1,3) = yAcorr(iMin,jMin);
busStatsMin(2,3) = cAcorr(iMin,jMin);
busStatsMin(3,3) = iAcorr(iMin,jMin);

griddz(iMin)
gridp11(jMin)
disp(busStatsMin)

