function [x] = bisection2(Z,x0)

TOL = 10e-6;

lb = 1.631317514299226; ub = 1.635; x = x0;
while abs(Z(x)) > TOL
    x = 0.5*(lb+ub);
    if Z(x) < 0
        ub = x;
    else
        lb = x;
    end;
    
    disp(x)
end;

end

