function [vfun1,vfun2,K1,K2,kSim,busStats] = everything(dz,p11,N,TOL,alpha,beta,delta,rho,kss,burn,T,gridk)

z1 = 1 + dz;
z2 = 1 - dz;
p12 = 1 - p11;
p21 = p12;
p22 = 1 - p21;

% initial guess of parameters of value function:

% state 1
dOld1 = 1;
gOld1 = 1;
hOld1 = 1;

% state 2
dOld2 = 1;
gOld2 = 1;
hOld2 = 1;

% state 1
dNew1 = 10;
gNew1 = 10;
hNew1 = 10;

% state 2
dNew2 = 10;
gNew2 = 10;
hNew2 = 10;

% These vectors hold the values of the vfuns at each grid point. 
% Used only for feeding to the polyfit function based on the parameter
% estimates.
vfun1 = zeros(1,N);
vfun2 = zeros(1,N);
K1 = zeros(1,N); % policy function in state z_1
K2 = zeros(1,N); % policy function in state z_2

% Quadratic Approximation of the value function using the nonlinear FOC to
% solve for k' at each point in the grid
it=1;
while abs(dNew1 - dOld1) + abs(gNew1 - gOld1) + abs(hNew1 - hOld1) + abs(dNew2 - dOld2) + abs(gNew2 - gOld2) + abs(hNew2 - hOld2) > TOL
    disp(it)
    
    for i=1:N
        % shock z_1
        k = gridk(i);
        Z = @(kp) (z1*k^alpha + (1-delta)*k - kp)^(-rho) - beta * (p11*(gOld1 + 2*hOld1*kp) + p12*(gOld2 + 2*hOld2*kp));
        kp = fzero(Z,k); % optimal value of k' to be plugged into value function
        K1(i) = kp; % store the value
        vfun1(i) = (1/(1-rho)) * (z1*k^alpha + (1-delta)*k - kp)^(1-rho) + beta * (p11*(dOld1 + gOld1*kp + hOld1*kp^2) + p12*(dOld2 + gOld2*kp + hOld2*kp^2));

        % shock z_2
        k = gridk(i);
        Z = @(kp) (z2*k^alpha + (1-delta)*k - kp)^(-rho) - beta * (p21*(gOld1 + 2*hOld1*kp) + p22*(gOld2 + 2*hOld2*kp));
        kp = fzero(Z,k);
        K2(i) = kp;
        vfun2(i) = (1/(1-rho)) * (z2*k^alpha + (1-delta)*k - kp)^(1-rho) + beta * (p21*(dOld1 + gOld1*kp + hOld1*kp^2) + p22*(dOld2 + gOld2*kp + hOld2*kp^2));
        
    end
    
    p1 = polyfit(gridk,vfun1,2); % fits polynomial of degree 2 to the data
    p2 = polyfit(gridk,vfun2,2); % fits polynomial of degree 2 to the data
    
    % Replace old parameter estimates with the current "new" ones
    dOld1 = dNew1;
    gOld1 = gNew1;
    hOld1 = hNew1;
    dOld2 = dNew2;
    gOld2 = gNew2;
    hOld2 = hNew2;
    
    % New parameter estimates
    dNew1 = p1(3);
    gNew1 = p1(2);
    hNew1 = p1(1);
    dNew2 = p2(3);
    gNew2 = p2(2);
    hNew2 = p2(1);

    it = it + 1;
end

% current state
k = kss; % start the simulation at the nonstochastic steady state
z = z1; % we start in state z_1

% holds final simulated values of capital in states 1 and 2, resp.
kSim = zeros(1,T);
ySim = zeros(1,T); y = 0;
cSim = zeros(1,T); c = 0;
iSim = zeros(1,T); i = 0;

for t=1:T
    % store the results of the simulation
    kSim(t) = k;
    ySim(t) = y;
    cSim(t) = c;
    iSim(t) = i;
    
    r = rand;
    if z==z1
        % determine the capital stock to transition to next
        Z = @(kp) (z1*k^alpha + (1-delta)*k - kp)^(-rho) - beta * (p11*(gNew1 + 2*hNew1*kp) + p12*(gNew2 + 2*hNew2*kp));
        tmp = k;
        k = max(fzero(Z,k), 10e-8); % the max ensures that capital is always positive
        y = z1*tmp^alpha;
        c = y + (1-delta)*tmp - k;
        i = y - c;

        % determine which state to transition to next:
        if r < p11
            z = z1;
        else
            z = z2;
        end
    else
        Z = @(kp) (z2*k^alpha + (1-delta)*k - kp)^(-rho) - beta * (p11*(gNew1 + 2*hNew1*kp) + p12*(gNew2 + 2*hNew2*kp));
        tmp = k;
        k = max(fzero(Z,k), 10e-8);
        y = z2*tmp^alpha;
        c = y + (1-delta)*tmp - k;
        i = y - c;

        if r < p21
            z = z1;
        else
            z = z2;
        end
    end
end

% Calculate the business cycle statistics when the economy is in the
% stochastic steady state.

busStats = zeros(3); % same format as the table in the homework

% Note that standard deviations are in percent
busStats(1,1) = 100 * std(ySim(burn+1:T));
busStats(2,1) = 100 * std(cSim(burn+1:T));
busStats(3,1) = 100 * std(iSim(burn+1:T));

busStats(1,2) = 1;
busStats(2,2) = corr(cSim(burn+1:T)',ySim(burn+1:T)');
busStats(3,2) = corr(iSim(burn+1:T)',ySim(burn+1:T)');

ySimLag = zeros(1,T-1);
cSimLag = zeros(1,T-1);
iSimLag = zeros(1,T-1);
for t=1:T-1
    ySimLag(t+1) = ySim(t);
    cSimLag(t+1) = cSim(t);
    iSimLag(t+1) = iSim(t);
end
busStats(1,3) = corr(ySim(burn+1:T)',ySimLag(burn+1:T)');
busStats(2,3) = corr(cSim(burn+1:T)',cSimLag(burn+1:T)');
busStats(3,3) = corr(iSim(burn+1:T)',iSimLag(burn+1:T)');

end

