clc;
clear all;
close all;

y = @(x) ((4-tan(x))/((3*(x^2))+sin(x)));
f = @(x)  (x^3 - 2*x*((4-tan(x))/((3*(x^2))+sin(x))) + ((4-tan(x))/((3*(x^2))+sin(x)))^7 - 4*(x^3)*((4-tan(x))/((3*(x^2))+sin(x))) - 5);

syms x;
F = (x^3 - 2*x*((4-tan(x))/((3*(x^2))+sin(x))) + ((4-tan(x))/((3*(x^2))+sin(x)))^7 - 4*(x^3)*((4-tan(x))/((3*(x^2))+sin(x))) - 5);
f_dash = diff(F);
f_dash = matlabFunction(f_dash);

i=0;
x = 1;
fprintf('n\tx(n)\t\t\tf(x(n))\n')
while (abs(f(x))>10e-08)
    X(i+1,1)=i;
    X(i+1,2)=x;
    X(i+1,3)=f(x);
    fprintf('%d\t%e\t\t%13e \t\n',i,x,f(x));
    x = x - f(x)/f_dash(x);
    i = i+1;
end
X(i+1,1)=i;
X(i+1,2)=x;
X(i+1,3)=f(x);
fprintf('%d\t%e\t\t%13e \t\n',i,x,f(x));