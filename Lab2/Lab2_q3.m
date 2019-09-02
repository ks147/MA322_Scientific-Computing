clear all;
close all;
clf;
clc;
f = @(x) x.^2 - 1;
FixedPoint(f,1,1,-1,1);
FixedPoint(f,2,2,-2,2);

f = @(x) 1 + 2.*x - x.^2;
FixedPoint(f,1,3,-1,1);
FixedPoint(f,2,4,-1,1);


f = @(x) 0.5.*(1 + 3.*x - x.^2);
FixedPoint(f,1,5,1,2);
FixedPoint(f,1,6,1,2);

function FixedPoint(f,x0,fig,a,b)

xx = a:0.0001:b;
yy = f(xx);
figure(fig);
plot(xx,yy,'Color','b','LineWidth',1);
hold on;
line([a,b],[a,b],'Color','g','LineWidth',1)
x=[];
x(1) = x0;
n=1;
fprintf('%s %20s %20s \n','n','x(n)','f(x(n))-x(n)');

while(abs(f(x(n))-x(n)) > 10^-3 && n<20)
    x(n+1) = f(x(n));
    fprintf('%d %20f %20f \n',n,x(n),f(x(n))-x(n));
    line([x(n) f(x(n))],[f(x(n)) f(x(n))],'Color','r','LineWidth',0.5);
    line([f(x(n)) x(n+1)],[f(x(n)) f(x(n+1))],'Color','r','LineWidth',0.5);
    n = n+1;
end
fprintf('%d %20f %20f \n\n',n,x(n),f(x(n))-x(n));
end
