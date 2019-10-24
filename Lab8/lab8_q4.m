clear all;
clc;
close all;

f = @(t,y) 1 + (y/t) + (y/t).^2;
alpha = 0;
actual_sol = @(t) t.*tan(log(t));
a = 1;
b = 3;
h = 0.002;
starting_value = actual_sol([a,a+h,a+2*h]);
Adams_Bashforth(a,b,h,f,starting_value,actual_sol,1);
Adams_Moulton(a,b,h,f,starting_value,actual_sol,2);

starting_value = RK_order4(a,a+2*h,h,alpha,f);
Adams_Bashforth(a,b,h,f,starting_value,actual_sol,3);
Adams_Moulton(a,b,h,f,starting_value,actual_sol,4);

f = @(t,y) sin(t) - y;
alpha = 1;
h = 0.005;
a = 0;
b = 1;
actual_sol = @(t) (sin(t) -cos(t))/2 +1.5.*exp(-t);
starting_value = actual_sol([a,a+h,a+2*h]);
Adams_Bashforth(a,b,h,f,starting_value,actual_sol,5);
Adams_Moulton(a,b,h,f,starting_value,actual_sol,6);

starting_value = RK_order4(a,a+2*h,h,alpha,f);
Adams_Bashforth(a,b,h,f,starting_value,actual_sol,7);
Adams_Moulton(a,b,h,f,starting_value,actual_sol,8);

function [t error] = Adams_Bashforth(a,b,h,f,starting_value,actual_sol,fig_no);
y = starting_value;
w = [23/12 -16/12 5/12];
t = [a:h:b];
n = length(t);

for i=3:n-1
    F = zeros(1,3);
    for j=1:3
        F(j) = f(t(i-j+1),y(i-j+1));
    end
    y(i+1) = y(i) + h*dot(w,F);
end
error = abs(y - actual_sol(t));
figure(fig_no)
plot(t,error);
title('Adams Bashforth');
xlabel('t')
ylabel('error');
end

function [t error] = Adams_Moulton(a,b,h,f,starting_value,actual_sol,fig_no)
y = starting_value;
w_bash = [23/12 -16/12 5/12];
w = [9/24 19/24 -5/24 1/24];

t = [a:h:b];
n = length(t);

for i=3:n-1
    F = zeros(1,3);
    for j=1:3
        F(j) = f(t(i-j+1),y(i-j+1));
    end
    y(i+1) = y(i) + h*dot(w_bash,F);  %Adam's Bashforth prediction
    F = zeros(1,4);
    for j=1:4
        F(j) = f(t(i-j+2),y(i-j+2));
    end
    y(i+1) = y(i) + h*dot(w,F); %Adam's Moulton correction
end
error = abs(y - actual_sol(t));
figure(fig_no)
plot(t,error);
title('Adams Moulton');
xlabel('t')
ylabel('error');
end

function y = RK_order4(a,b,h,alpha,f)
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

end