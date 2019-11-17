clc;
clf;
close all;
clear all;
%Q5(a)

fprintf('Q5(a) \n');
Newton(1,2,10^-5,'a');

%Q5(b)
fprintf('Q5(b) \n');
Newton(1.3,2,10^-5,'b')

%Q5(c)
fprintf('Q5(c) \n');
Newton(2,3,10^-5,'c');
Newton(3,4,10^-5,'c');

%Q5(d)
fprintf('Q5(d) \n');
Newton(1,2,10^-5,'d');
Newton(exp(1),4,10^-5,'d');

%Q5(e)
fprintf('Q5(e) \n');
Newton(0,1,10^-5,'e');
Newton(3,5,10^-5,'e');

%Q5(f)
fprintf('Q5(f) \n');
Newton(0,1,10^-5,'f');
Newton(3,4,10^-5,'f');
Newton(6,7,10^-5,'f');



function y = f_(x,part)
   if(part=='a')
        f = @(x) exp(x) + 2^-x + 2*cos(x) - 6;
    elseif(part=='b')
        f = @(x) log(x-1) + cos(x-1);
    elseif(part=='c')
        f = @(x) 2*x*cos(2*x)-(x-2)^2;
    elseif(part=='d')
        f = @(x) (x-2)^2 - log(x);
    elseif(part=='e')
        f = @(x) exp(x) - 3*(x^2);
    else
        f = @(x) sin(x) - exp(-x);
    end
    h = 10^-7;
    y = (f(x+h)-f(x-h))/(2*h);
end
function Newton(a,b,limit,part)
    if(part=='a')
        f = @(x) exp(x) + 2^-x + 2*cos(x) - 6
    elseif(part=='b')
        f = @(x) log(x-1) + cos(x-1)
    elseif(part=='c')
        f = @(x) 2*x*cos(2*x)-(x-2)^2
    elseif(part=='d')
        f = @(x) (x-2)^2 - log(x)
    elseif(part=='e')
        f = @(x) exp(x) - 3*(x^2)
    else
        f = @(x) sin(x) - exp(-x)
    end
            
    x=[];
    x(1) = (a+b)/2;
    i=1;
    fprintf('%f <= x <= %f \n',a,b);
    fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');

    while(abs(f(x(i))) > limit)
        fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
        x(i+1) = x(i) - (f(x(i))/f_(x(i),part));
        i=i+1;
    end
    fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
end
