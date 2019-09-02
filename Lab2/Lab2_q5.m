clear all;
clf;
clc;
close all;
%Q5(a)
fprintf('Q5(a)\nx^3 - x - 2\n');
Muller(1,1.2,1.4,'a');
%Q5(b)
fprintf('Q5(b)\n1 + 2x - tanx\n');
Muller(1.5,1.4,1.3,'b');

function y = f2(x1,x2,part)
     if(part=='a')
        f = @(x) x^3 - x - 2;
    else
        f = @(x) 1 + 2*x  - tan(x);
    end
    y = (f(x2)-f(x1))/(x2-x1);
end
function y = f3(x1,x2,x3,part)
    y = (f2(x2,x3,part)-f2(x1,x2,part))/(x3-x1);     
end

function Muller(a,b,c,part)
    if(part=='a')
        f = @(x) x^3 - x - 2;
    else
        f = @(x) 1 + 2*x  - tan(x);
    end
x = [];
x(1) = a;
x(2) = b;
x(3) = c;
i = 3;
fprintf('%s %20s %20s \n','n','x(n)','f(x(n))');
while(abs(f(x(i))) > 10^-5)
    w = f2(x(i),x(i-1),part) + f2(x(i),x(i-2),part) - f2(x(i-1),x(i-2),part);
    D = w^2 - 4*f(x(i))*f3(x(i),x(i-1),x(i-2),part);
    if(abs(w+sqrt(D)) > abs(w-sqrt(D)))
        den = w+sqrt(D);
    else
        den = w - sqrt(D);
    end
    x(i+1) = x(i) - (2*f(x(i))/den);
    fprintf('%d %20f %20f \n',i,x(i),f(x(i)));
    i = i+1;
end
fprintf('%d %20f %20f \n\n',i,x(i),f(x(i)));
end