clear all;
clf;
clc;

%(a)
lambda = -5;
f = @(t,y) lambda*y + (1-lambda)*cos(t) - (1+lambda)*sin(t);
actual_y = @(t) sin(t)+cos(t);
a = 0;
b = 10;
h = 1e-2;
alpha = 1;
fprintf('lambda = -5, h = 0.01\n');
explicit_euler(a,b,f,h,alpha,actual_y,1);
implicit_euler(a,b,f,h,alpha,actual_y,2);
%lambda = -5 ,h = 0.001
h = 1e-3;
fprintf('lambda = -5, h= 0.001\n');
explicit_euler(a,b,f,h,alpha,actual_y,3);
implicit_euler(a,b,f,h,alpha,actual_y,4);

%lambda = -5 h = 0.005
h = 5*1e-3;
fprintf('lambda = -5, h= 0.005\n');
explicit_euler(a,b,f,h,alpha,actual_y,5);
implicit_euler(a,b,f,h,alpha,actual_y,6);
%lambda = -5 h = 0.0025
h = 0.0025;
fprintf('lambda = -5, h= 0.0025\n');
explicit_euler(a,b,f,h,alpha,actual_y,7);
implicit_euler(a,b,f,h,alpha,actual_y,8);

%(b)
%lambda = 5 ,h = 0.01
h = 0.01;
lambda = 5;
f = @(t,y) lambda*y + (1-lambda)*cos(t) - (1+lambda)*sin(t);
fprintf('lambda = 5, h= 0.01\n');
explicit_euler(a,b,f,h,alpha,actual_y,9);
implicit_euler(a,b,f,h,alpha,actual_y,10);

%lamda = 5 , h = 0.001
h = 0.001;
fprintf('lambda = 5, h= 0.001\n');
explicit_euler(a,b,f,h,alpha,actual_y,11);
implicit_euler(a,b,f,h,alpha,actual_y,12);

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

