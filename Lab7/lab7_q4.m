clear all;
clf;
clc;
close all;
%(a)
fprintf('Q4(a)\n');
f = @(t,y) (2-2*t*y)/(t.^2+1);
a = 0;
b = 1;
alpha = 1;
h = 0.01;
actual_y = @(t) (2*t-1)/(t.^2+1);
explicit_euler(a,b,f,h,alpha,actual_y,1);
implicit_euler(a,b,f,h,alpha,actual_y,2);
trapezoid_method(a,b,f,h,alpha,actual_y,3);
modified_euler(a,b,f,h,alpha,actual_y,4);

%(b)
fprintf('Q4(b)\n');
f = @(t,y) (y.^2 + y)/t;
a = 1;
b = 3;
alpha = -2;
h = 0.02;
actual_y = @(t) 2*t/(1-2*t);
explicit_euler(a,b,f,h,alpha,actual_y,5);
implicit_euler(a,b,f,h,alpha,actual_y,6);
trapezoid_method(a,b,f,h,alpha,actual_y,7);
modified_euler(a,b,f,h,alpha,actual_y,8);

%(c)
fprintf('Q4(c)\n');
f = @(t,y) 1 + (y/t) + (y/t).^2;
a = 1;
b = 3;
alpha = 0;
h = 0.002;
actual_y = @(t) t.*tan(log(t));
explicit_euler(a,b,f,h,alpha,actual_y,9);
implicit_euler(a,b,f,h,alpha,actual_y,10);
trapezoid_method(a,b,f,h,alpha,actual_y,11);
modified_euler(a,b,f,h,alpha,actual_y,12);

%(d)
fprintf('Q4(d)\n');
f = @(t,y) exp(t-y);
a = 0;
b = 1;
alpha = 1;
h = 0.005;
actual_y = @(t) log(exp(t) + exp(1) -1);
explicit_euler(a,b,f,h,alpha,actual_y,13);
implicit_euler(a,b,f,h,alpha,actual_y,14);
trapezoid_method(a,b,f,h,alpha,actual_y,15);
modified_euler(a,b,f,h,alpha,actual_y,16);


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

function [y,error,t] = implicit_euler(a,b,f,h,alpha,actual_y,fig_no)
t = [a:h:b];
y(1) = alpha;
N = length(t);
error(1) = abs(y(1)-actual_y(t(1)));
%fprintf('%s\t\t%s\t%s\t\t%s\n','t(i)','y_appx(i)','y_actual(t(i))','error(i)');
%fprintf('%f\t%f\t%f\t%f\n',t(1),y(1),actual_y(t(1)),error(1));
for i=2:N
    y_prev = y(i-1);
    y_new = y_prev + h*f(t(i),y_prev);
        
    while abs(y_new-y_prev) > 1e-8
        save = y_new;
        y_new = y_prev + h*f(t(i),y_prev);
        y_prev =  save;
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

function [y,error,t] = trapezoid_method(a,b,f,h,alpha,actual_y,fig_no)

t = [a:h:b];
y(1) = alpha;
N = length(t);
error(1) = abs(y(1)-actual_y(t(1)));
%fprintf('%s\t\t%s\t%s\t\t%s\n','t(i)','y_appx(i)','y_actual(t(i))','error(i)');
%fprintf('%f\t%f\t%f\t%f\n',t(1),y(1),actual_y(t(1)),error(1));
for i=2:N
    y_exp = y(i-1)+h*f(t(i-1),y(i-1));
    y(i) = y(i-1) + 0.5*h*(f(t(i-1),y(i-1)) + f(t(i),y_exp));
    error(i) = abs(y(i)-actual_y(t(i)));
    %fprintf('%f\t%f\t%f\t%f\n',t(i),y(i),actual_y(t(i)),error(i));
end
figure(fig_no);
subplot(1,2,1);
plot(t,error);
xlabel('t(i)');
ylabel('error(i)');
title('Trapezoid Method');

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