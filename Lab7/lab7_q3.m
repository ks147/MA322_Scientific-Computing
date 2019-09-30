clear all;
clc;
clf;
close all;

lambda = -1;
f = @(t,y) lambda*y + 1/(1+t^2) - lambda*atan(t);
actual_y = @(t) atan(t);
alpha = 0;
a = 0;
b = 1;
%lambda = -1,h = 0.01
h = 0.01;
explicit_euler(a,b,f,h,alpha,actual_y,1);
implicit_euler(a,b,f,h,alpha,actual_y,2);
modified_euler(a,b,f,h,alpha,actual_y,3);

%lambda = -1,h = 0.001
h = 0.001;
explicit_euler(a,b,f,h,alpha,actual_y,4);
implicit_euler(a,b,f,h,alpha,actual_y,5);
modified_euler(a,b,f,h,alpha,actual_y,6);

%lambda = -10,h = 0.01
lambda = -10;
f = @(t,y) lambda*y + 1/(1+t^2) - lambda*atan(t);
h = 0.01;
explicit_euler(a,b,f,h,alpha,actual_y,7);
implicit_euler(a,b,f,h,alpha,actual_y,8);
modified_euler(a,b,f,h,alpha,actual_y,9);

%lambda = -10, h = 0.001
h = 0.001;
explicit_euler(a,b,f,h,alpha,actual_y,10);
implicit_euler(a,b,f,h,alpha,actual_y,11);
modified_euler(a,b,f,h,alpha,actual_y,12);

%lambda = -50,h = 0.01
lambda = -50;
f = @(t,y) lambda*y + 1/(1+t^2) - lambda*atan(t);
h = 0.01;
explicit_euler(a,b,f,h,alpha,actual_y,13);
implicit_euler(a,b,f,h,alpha,actual_y,14);
modified_euler(a,b,f,h,alpha,actual_y,15);

%lambda = -50
h = 0.001;
explicit_euler(a,b,f,h,alpha,actual_y,16);
implicit_euler(a,b,f,h,alpha,actual_y,17);
modified_euler(a,b,f,h,alpha,actual_y,18);


%%functions 
function [y,error,t] = explicit_euler(a,b,f,h,alpha,actual_y,fig_no)
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
figure(fig_no);
subplot(1,2,1);
plot(t,error);
xlabel('t(i)');
ylabel('error(i)');
title('Explicit Euler');

subplot(1,2,2);
plot(t,y,'--','linewidth',2);
hold on;
plot(t,actual_y(t),'linewidth',2);
hold off;
xlabel('t');
ylabel('y(t)');
title('Appx. Y(t) vs Actual Y(t)');

fprintf('\n');
end

function implicit_euler(a,b,f,h,alpha,actual_y,fig_no)
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
        y_new = y_prev + h*f(t(i-1),y_prev);
        diff = abs(y_new - y_prev);
        y_prev =  y_new;
    end
    y(i) = y_new;
    error(i) = abs(y(i)-actual_y(t(i)));
    %fprintf('%f\t%f\t%f\t%f\n',t(i),y(i),actual_y(t(i)),error(i));
end
figure(fig_no);
subplot(1,2,1);
plot(t,error);
xlabel('t(i)');
ylabel('error(i)');
title('Implicit Euler');

subplot(1,2,2);
plot(t,y,'--','linewidth',2);
hold on;
plot(t,actual_y(t),'linewidth',2);
hold off;
xlabel('t');
ylabel('y(t)');
title('Appx. Y(t) vs Actual Y(t)');
fprintf('\n');
end

function [y,error,t] = modified_euler(a,b,f,h,alpha,actual_y,fig_no)

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
figure(fig_no);
subplot(1,2,1);
plot(t,error);
xlabel('t(i)');
ylabel('error(i)');
title('Modified Euler');

subplot(1,2,2);
plot(t,y,'--','linewidth',2);
hold on;
plot(t,actual_y(t),'linewidth',2);
hold off;
xlabel('t');
ylabel('y(t)');
title('Appx. Y(t) vs Actual Y(t)');
fprintf('\n');
end