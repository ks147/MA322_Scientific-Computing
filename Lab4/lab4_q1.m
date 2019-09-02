clear all;
clc;
close all;


f = @(x) exp(x);

fprintf('Forward Diff Interpolation\n');
Forward_Diff_Interpolation(f,10,0.5,1,2.25);
fprintf('Backward Diff Interpolation\n');
Backward_Diff_Interpolation(f,30,0.05,1,2.25);
function Forward_Diff_Interpolation(f,n,h,x0,eval_at)
    
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
        p = p + (u_cal_for(mu,i-1)*forward_diff(1,i))/factorial(i-1);
    end;
    P = matlabFunction(p);
    Error = @(x) (abs(P(x)-f(x)));
    fprintf('value at:%e = %e\nerror:%e\n\n',eval_at,P(eval_at),Error(eval_at));
    XX = [x(1):1e-4:x(n)];
    figure(1);
    plot(XX,Error(XX));
end

function y = u_cal_for(u,n)
    y = u;
    for i=1:n-1
        y = y*(u-i);
    end
end
function Backward_Diff_Interpolation(f,n,h,x0,eval_at)
    back_diff = zeros(n,n);
    x = [];
    x(1) = x0;
    
    for i=2:n
        x(i) = x(i-1) + h; 
    end;
    
    for i=1:n
        back_diff(i,1) = f(x(i));
    end;

    for i= 1:n-1
        j = n-1;
        while j >= i
            back_diff(j+1,i+1) = back_diff(j+1,i) - back_diff(j,i);
            j = j-1;
        end
    end
    syms z;
    %Evaluate the polynomial at x = 2.25
    mu = (z-x(n))/h;

    p = back_diff(n,1);
    
    for i = 2:n
        p = p + (u_cal_back(mu,i-1)*back_diff(n,i))/factorial(i-1);
    end;
    P = matlabFunction(p);
    Error = @(x) (abs(P(x)-f(x)));
    fprintf('value at:%e = %e\nerror:%e\n\n',eval_at,P(eval_at),Error(eval_at));
    XX = [x(1):1e-4:x(n)];
    figure(2);
    plot(XX,Error(XX));
end
function y = u_cal_back(u,n)
    y = u;
    for i=1:n-1
        y = y*(u+i);
    end
end
