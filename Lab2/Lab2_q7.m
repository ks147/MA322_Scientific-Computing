clear all;
clf;
clc;
close all;
%Q7(a)
fprintf('Roots of z^4 - 2z^3 - 2iz^2 + 4iz = 0\n');
Muller(0,0.1,0.2,'a');
Muller(-1,-1.2,-1.3,'a');
Muller(1.4,1.5,1.6,'a');
Muller(3.4,3.5,3.6,'a');

%Q7(b)
fprintf('Roots of z = e^z\n');
Muller(0.7,0.75,0.8,'b');

function y = f2(x1,x2,part)
     if(part=='a')
        f = @(z) z^4 - 2*z^3 - 2i*z^2 + 4i*z;
    else
        f = @(z) exp(z) - z;
    end
    y = (f(x2)-f(x1))/(x2-x1);
end
function y = f3(x1,x2,x3,part)
    y = (f2(x2,x3,part)-f2(x1,x2,part))/(x3-x1);     
end

function Muller(a,b,c,part)
     if(part=='a')
         f = @(z) z^4 - 2*z^3 - 2i*z^2 + 4i*z;
    else
        f = @(x) x^2 - sin(x);
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
    fprintf('%d %20f+%fj %20f+%fj \n',i,real(x(i)),imag(x(i)),real(f(x(i))),imag(f(x(i))));
    
    i = i+1;
end
fprintf('%d %20f+%fj %20f+%fj \n\n',i,real(x(i)),imag(x(i)),real(f(x(i))),imag(f(x(i))));
end