clear all;
clc;
close all;

f = @(x) (cos(x+((2)^0.5)) + x*(x/2 + ((2)^0.5)));

syms x
F = (cos(x+((2)^0.5)) + x*(x/2 + ((2)^0.5)));
X1 = Newton(f,F,-1);
p = multiplicity(F,-sqrt(2));
Y1 = modifNewton(f,F,-1,p);

f = @(x)(exp(6*x)+3*((log(2))^2)*exp(2*x)-(log(8)*exp(4*x))-(log(2))^3);
syms x
F = (exp(6*x)+3*((log(2))^2)*exp(2*x)-(log(8)*exp(4*x))-(log(2))^3);
X2 = Newton(f,F,0);
p = multiplicity(F,-1.833);
Y2 = modifNewton(f,F,0,p);

function p = multiplicity(F,root)
    F = diff(F);
    df = matlabFunction(F);
    p = 1;
    while(abs(df(root))<1e-2)
        F = diff(F);
        df = matlabFunction(F);
        p = p+1;
    end
end

function  X = Newton(f,F,x)
    F = diff(F);
    i = 0;
    df = matlabFunction(F);
    fprintf('\nn\tx(n)\t\t\tf(x(n))\n\n')
    while(abs(f(x))>1e-15)
        X(i+1,1)=i;
        X(i+1,2)=x;
        X(i+1,3)=f(x);
        fprintf('%d\t%e\t\t%13e \t\n',i,x,f(x));
        x = x-(f(x)/df(x));
        i = i+1;
    end
    X(i+1,1)=i;
    X(i+1,2)=x;
    X(i+1,3)=f(x);
    fprintf('%d\t%e\t\t%13e \t\n',i,x,f(x));
end

function X = modifNewton(f,F,x,p)
    F = diff(F);
    i = 0;
    df = matlabFunction(F);
    fprintf('\nn\tx(n)\t\t\tf(x(n))\n\n')
    while(abs(f(x))>1e-15)
        X(i+1,1)=i;
        X(i+1,2)=x;
        X(i+1,3)=f(x);
        fprintf('%d\t%e\t\t%13e \t\n',i,x,f(x));
        x = x-p*(f(x)/df(x));
        i = i+1;
    end
    X(i+1,1)=i;
    X(i+1,2)=x;
    X(i+1,3)=f(x);
    fprintf('%d\t%e\t\t%13e \t\n',i,x,f(x));
end