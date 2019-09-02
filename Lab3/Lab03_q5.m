clear all;
clc;
close all;

f = @(x) ((x^5)-2*(x^4)-2*(x^3)+8*(x^2)-7*x+2);
syms x
F = ((x^5)-2*(x^4)-2*(x^3)+8*(x^2)-7*x+2);
X1 = modifNewton(f,F,1.3,4,1);

f = @(x) ((x^5)-8*(x^4)+25*(x^3)-38*(x^2)+28*x-8);
syms x
F = ((x^5)-8*(x^4)+25*(x^3)-38*(x^2)+28*x-8);
X21 = modifNewton(f,F,0,2,1);
X22 = modifNewton(f,F,3,3,2);

function X=modifNewton(f,F,x,p,root)
    F = diff(F);
    i = 0;
    df = matlabFunction(F);
    fprintf('\nn\tx(n)\t\tf(x(n))\t\te(n)\t\te(n+1)\t\te(n+1)/e(n)\tlog(e(n+1)/e(n))\t\t\t\n\n')
    while(abs(f(x))>1e-6)
        X(i+1,1)=i;
        X(i+1,2)=x;
        X(i+1,3)=f(x);
        X(i+1,4)=abs(x-root);
        xx = x-p*(f(x)/df(x));
        X(i+1,5)=abs(xx-root);
        X(i+1,6)=abs(xx-root)/abs(x-root);
        X(i+1,7)=log(abs(xx-root)/abs(x-root));
        fprintf('%d\t%e\t%13e\t%13e\t%13e\t%13e\t%13e\n',i,x,f(x),abs(x-root),abs(xx-root),abs(xx-root)/abs(x-root),log(abs(xx-root)/abs(x-root)));
        x = x-p*(f(x)/df(x));
        i = i+1;
    end
    X(i+1,1)=i;
    X(i+1,2)=x;
    X(i+1,3)=f(x);
    X(i+1,4)=abs(x-root);
    if(df(x)==0 && f(x)==0 )
        xx=x
    end
    X(i+1,5)=abs(xx-root);
    X(i+1,6)=abs(xx-root)/abs(x-root);
    X(i+1,7)=log(abs(xx-root)/abs(x-root));
    fprintf('%d\t%e\t%13e\t%13e\t%13e\t%13e\t%13e\n',i,x,f(x),abs(x-root),abs(xx-root),abs(xx-root)/abs(x-root),log(abs(xx-root)/abs(x-root)));
end