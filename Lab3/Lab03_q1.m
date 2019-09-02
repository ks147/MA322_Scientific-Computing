clf;
clc;
clear all;
close all;
%Q1(a)
f1 = @(x,y) x*(y.^2) + x^2*y + x^4 - 3;
f2 = @(x,y) (x^3)*y^5 - 2*(x^5*y) - x^2 +2;
a = [1;1];
Newton(a,f1,f2);
f1 = @(x,y) sin(4*pi*x*y)-2*y-x;
f2 = @(x,y) (4*pi-1)*(exp(2*x)-exp(1))/(4*pi) + 4*exp(1)*y^2 -2*exp(1)*x ;
a = [0;0];
Newton(a,f1,f2);
%
function Newton(a,f1,f2)
    syms x y
    F = jacobian([f1(x,y),f2(x,y)],[x y]);
    fprintf('%s %20s %20s\n','n','x(1)','x(2)');   
    i = 1;
    for i=1:3
        fprintf('%d %20f %20f\n',i,a(1),a(2));
        A = subs(F,{x,y},{a(1),a(2)});
        f = zeros(2,1);
        f(1) = f1(a(1),a(2));
        f(2) = f2(a(1),a(2));
        delta = -A\f;
        a = a+delta;
    end
    fprintf('%d %20f %20f\n',i+1,a(1),a(2));
end