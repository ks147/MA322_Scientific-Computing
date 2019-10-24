clear all;
clc;

%Runge Kutta order 4

f = @(t,y) F(t,y);

alpha = 0;
a = 0;
b = 1;
h = 2^-7;
f(0.5*h,0)
actual_sol = @(t) t*asin(t);
RK_order4(a,b,h,alpha,f,actual_sol);

function RK_order4(a,b,h,alpha,f,actual_sol)
y = [];
y(1) = alpha;
w = [1/6,1/3,1/3,1/6];
t = [a:h:b];
n = length(t);
c = [0,1/2,1/2,1];

error(1) = 0;
for i=1:n-1
    k = zeros(1,4);
    k(1) = h*f(t(i),y(i));
    k(2) = h*f(t(i)+c(2)*h,y(i)+0.5*k(1));
    k(3) = h*f(t(i)+c(3)*h,y(i)+0.5*k(2));
    k(4) = h*f(t(i)+c(4)*h,y(i)+k(3));
    y(i+1) = y(i) + dot(w,k);
    error(i+1) = abs(y(i+1) - actual_sol(t(i+1)));
end
%plot solution
plot(t,error);
end

function x = F(t,y)

if t==0 && y==0
x = 0;
else
x = (y/t) + t*sec(y/t);
end
end



