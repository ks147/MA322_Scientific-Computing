clear all;
clc;
clf;
close all;

f = @(t,y) t+y;
actual_y = @(t) 2*exp(t) - t - 1;
alpha = 1;
a = 0;
b = 1;

% h = 0.05
h = 0.05;
[y1,error1,t1] = explicit_euler(a,b,f,h,alpha,actual_y);
[y2,error2,t2] = modified_euler(a,b,f,h,alpha,actual_y);
n = length(t1);
fprintf('h = %f\n',h);
fprintf('Euler method gives y(1) = %f\n',y1(n));
fprintf('Modified Euler method gives y(1) = %f\n',y2(n));
fprintf('Actual solution gives y(1) = %f\n\n',actual_y(1));

figure(1);
subplot(1,2,1);
plot(t1,error1);
xlabel('t');
ylabel('Error');
title('Euler method');

subplot(1,2,2);
plot(t2,error2);
xlabel('t');
ylabel('Error');
title('Modified Euler method');

figure(2);
subplot(1,2,1);
hold on;
plot(t1,y1);
plot(t1,actual_y(t1));
xlabel('t');
ylabel('y(t)');
title('Euler method Y(t) vs Actual Y(t)');

subplot(1,2,2);
hold on;
plot(t2,y2);
plot(t2,actual_y(t2));
xlabel('t');
ylabel('y(t)');
title('Modified Euler method Y(t) vs Actual Y(t)');

%h = 0.005
h = 0.005;
[y1,error1,t1] = explicit_euler(a,b,f,h,alpha,actual_y);
[y2,error2,t2] = modified_euler(a,b,f,h,alpha,actual_y);
n = length(t1);
fprintf('h = %f\n',h);
fprintf('Euler method gives y(1) = %f\n',y1(n));
fprintf('Modified Euler method gives y(1) = %f\n',y2(n));
fprintf('Actual solution gives y(1) = %f\n',actual_y(1));

figure(3);
subplot(1,2,1);
plot(t1,error1);
xlabel('t');
ylabel('Error');
title('Euler method');

subplot(1,2,2);
plot(t2,error2);
xlabel('t');
ylabel('Error');
title('Modified Euler method');

figure(4);
subplot(1,2,1);
hold on;
plot(t1,y1);
plot(t1,actual_y(t1));
xlabel('t');
ylabel('y(t)');
title('Euler method Y(t) vs Actual Y(t)');

subplot(1,2,2);
hold on;
plot(t2,y2);
plot(t2,actual_y(t2));
xlabel('t');
ylabel('y(t)');
title('Modified Euler method Y(t) vs Actual Y(t)');


%%functions
function [y,error,t] = explicit_euler(a,b,f,h,alpha,actual_y)
t = [a:h:b];
y(1) = alpha;
N = length(t);
error(1) = abs(y(1)-actual_y(t(1)));
%fprintf('%s\t\t%s\t%s\t\t%s\n','t(i)','y_appx(i)','y_actual(t(i))','error(i)');
%fprintf('%f\t%f\t%f\t%f\n',t(1),y(1),actual_y(t(1)),error(1));
for i=2:N
    y(i) = y(i-1) + h*f(t(i-1),y(i-1));
    error(i) = abs(y(i)-actual_y(t(i)));
    %fprintf('%f\t%f\t%f\t%f\n',t(i),y(i),actual_y(t(i)),error(i));
end

end

function [y,error,t] = modified_euler(a,b,f,h,alpha,actual_y)

t = [a:h:b];
y(1) = alpha;
N = length(t);
error(1) = abs(y(1)-actual_y(t(1)));
%fprintf('%s\t\t%s\t%s\t\t%s\n','t(i)','y_appx(i)','y_actual(t(i))','error(i)');
%fprintf('%f\t%f\t%f\t%f\n',t(1),y(1),actual_y(t(1)),error(1));
for i=2:N
    y_prev = y(i-1);
    y_new = y(i-1)+10;
    diff = abs(y_new-y_prev);
    
    while(diff > 1e-8)
        y_new = y(i-1) + 0.5*h*(f(t(i-1),y(i-1)) + f(t(i),y_prev));
        diff = abs(y_new - y_prev);
        y_prev =  y_new;
    end
    y(i) = y_new;
    error(i) = abs(y(i)-actual_y(t(i)));
    %fprintf('%f\t%f\t%f\t%f\n',t(i),y(i),actual_y(t(i)),error(i));
end

end