clc;
clear all;
kappa = 6.22*10^-19;
n1 = 2*10^3;
n2 = 2*10^3;
n3 = 3*10^3;
f = @(t,x) kappa*((n1 - x/2).^2)*((n2 - x/2).^2)*(n3 - 3*x/4).^3;

alpha = 0;
a = 0;
b = 0.2;
h = 10^-4;

x = RK_order4(a,b,h,alpha,f);
fprintf('Amount of potassium hydroxide after 0.2s = %e\n',x);

function x = RK_order4(a,b,h,alpha,f)
y = [];
y(1) = alpha;
w = [1/6,1/3,1/3,1/6];
t = [a:h:b];
n = length(t);
c = [0,1/2,1/2,1];
k = [];

for i=1:n-1
    k = zeros(1,4);
    k(1) = h*f(t(i),y(i));
    k(2) = h*f(t(i)+c(2)*h,y(i)+0.5*k(1));
    k(3) = h*f(t(i)+c(3)*h,y(i)+0.5*k(2));
    k(4) = h*f(t(i)+c(4)*h,y(i)+k(3));
    y(i+1) = y(i) + dot(w,k);
end
x = y(n);
end