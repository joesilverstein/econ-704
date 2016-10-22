function [h] = bisection(Z)

TOL = 10e-6;

if Z(0) < 0 
    h = 0;
    return;
elseif Z(1) > 0
    h = 1;
    return;
else
    lb = 10e-8; ub = 1; h = 0.5;
    while abs(Z(h)) > TOL
        h = 0.5*(lb+ub);
        if Z(h) < 0
            ub = h;
        else
            lb = h;
        end;
    end;
end;

end

