close all;
clc;

f = @(x) atan(x);
e = [];
e(1) = 1;
for i = 2:33
    e(i) = e(i-1)+(8-1)/33;
end

Forward_Diff_Interpolation(f,11,5/11,1,e);
function Forward_Diff_Interpolation(f,n,h,x0,e)
    
    forward_diff = zeros(n,n);
    x = [];
    x(1) = x0;
    
    for i=2:n
        x(i) = x(i-1) + h; 
    end;
    for i=1:n
        forward_diff(i,1) = f(x(i));
    end;

    for i= 2:n
        for j=1:n-i+1
            forward_diff(j,i) = forward_diff(j+1,i-1) - forward_diff(j,i-1);
        end
    end
    
    syms z;
    mu = (z - x(1))/h;

    p = f(x(1));
    
    for i = 2:n
        p = p + (u_cal_for(mu,i-1)*(forward_diff(1,i))/factorial(i-1));
    end;
    fprintf('Newton Polynomial:%s\n',char(p));
    P = matlabFunction(p);
    Error = @(x) (abs(P(x)-f(x)));
    fprintf('n\tx(i)\t\tError(x(i))\n')
    for i=1:length(e)
        fprintf('%d\t%f\t%e\n',i,e(i),Error(e(i)));
    end
end

function y = u_cal_for(u,n)
    y = u;
    for i=1:n-1
        y = y*(u-i);
    end
end