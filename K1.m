function [k1] = K1(kT1,kT,T,E)
% K1 Given the last two values of the transition path, returns the first

kt = kT;
kt1 = kT1;
for i=1:T-2
    Z = @(kt2) E(kt2,kt1,kt);
    kt2 = fzero(Z,kt1); % solve for kT-2, using kT-1 as the initial guess
    kt = kt1;
    kt1 = max(kt2,10e-4);
end

k1 = kt2;

end

