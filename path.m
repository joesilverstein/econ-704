function [kPath] = path(kT1,kT,T,E)
% Given last two values of the transition path, returns the paths of
% k,h,c,y.

kPath = zeros(T,1);

kt = kT;
kt1 = kT1;
kPath(T) = kT;
kPath(T-1) = kT1;
for i=1:T-2
    t = T - i - 1; % backwards
    Z = @(kt2) E(kt2,kt1,kt);
    kt2 = fzero(Z,kt1); % solve for kT-2, using kT-1 as the initial guess
    kPath(t) = kt2;
    kt = kt1;
    kt1 = kt2;
end

end

