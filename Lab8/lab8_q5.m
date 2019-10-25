clear all;
clc;
close all;

f = @(t,y) -50.*y + 51.*cos(t) + 49.*sin(t);
alpha = 1;
a = 0;
b = 10;
actual_sol = @(t) sin(t) + cos(t);
h = 0.01;
starting_value = actual_sol([a,a+h,a+2*h]);
Adams_Bashforth(a,b,h,f,starting_value,actual_sol,1);
Adams_Moulton(a,b,h,f,starting_value,actual_sol,2);

h = 0.001;
starting_value = actual_sol([a,a+h,a+2*h]);
Adams_Bashforth(a,b,h,f,starting_value,actual_sol,3);
Adams_Moulton(a,b,h,f,starting_value,actual_sol,4);



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

figure(fig_no);
subplot(1,2,1);
plot(t,error);
title('Adams Bashforth');
xlabel('t')
ylabel('error');

subplot(1,2,2);
plot(t,y);
title('Adam Bashforth');
xlabel('t');
ylabel('y(t)');
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
figure(fig_no);
subplot(1,2,1);
plot(t,error);
title('Adams Moulton');
xlabel('t')
ylabel('error');

subplot(1,2,2);
plot(t,y);
title('Adam Moulton');
xlabel('t');
ylabel('y(t)');

end
