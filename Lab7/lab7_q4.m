clear all;
clf;
clc;
close all;

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