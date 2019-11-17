clc;
clf;
close all;
clear all;
%Q7(a)
fprintf('\n Newtons Method \n');
fprintf('first root \n');
Newton(-1,0,10^-6);
fprintf('\n Second root \n');
Newton(0,1.2,10^-6);

%Q7(b)
fprintf('\n Secant Method \n');
fprintf('first root \n');
Secant(-1,0.1,10^-6);
fprintf('\n Second root \n');
Secant(0,1,10^-6);


function y = f_(x)
    
    f = @(x) 230*(x^4) + 18*(x^3) + 9*(x^2) - 221*x - 9;
    h = 10^-7;
    y = (f(x+h)-f(x-h))/(2*h);
end
function Newton(a,b,limit)

    f = @(x) 230*(x^4) + 18*(x^3) + 9*(x^2) - 221*x - 9;
    x=[];
    x(1) = (a+b)/2;
    i=1;
    fprintf('%f <= x <= %f \n',a,b);
    fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');

    while(abs(f(x(i))) > limit)
        fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
        x(i+1) = x(i) - (f(x(i))/f_(x(i)));
        i=i+1;
    end
    fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
end

function Secant(a,b,limit)

    f = @(x) 230*(x^4) + 18*(x^3) + 9*(x^2) - 221*x - 9;
    x=[];
    x(1) = b;
    x(2) = a;
    i=2;
    fprintf('%f <= x <= %f \n',a,b);
    fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');

    while(abs(f(x(i))) > limit)
        fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
        derivative = f(x(i))-f(x(i-1));
        derivative = derivative/(x(i)-x(i-1));
        x(i+1) = x(i) - (f(x(i))/derivative);
        i=i+1;
    end
    fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
end

