clc;
clear all;

f = @(t,y) -y + (t.^(0.1)).*(1.1 + t);
actual_sol = @(t) t.^(1.1);
a = 0;
b = 5;
alpha = 0;
h = 0.005;
H = 1;
%Q1(a)
fprintf('\nQ1(a)\n');
RK_order2(a,b,h,H,alpha,f,actual_sol);
fprintf('\n');
%Q1(b)
h = 0.025;
fprintf('\nQ1(b)\n');
RK_order2(a,b,h,H,alpha,f,actual_sol);
%Q1(c)
h = 0.0125;
fprintf('\nQ1(c)\n');
RK_order2(a,b,h,H,alpha,f,actual_sol);

%Q1(d)
h = 0.00625;
fprintf('\nQ1(d)\n');
RK_order2(a,b,h,H,alpha,f,actual_sol);


function RK_order2(a,b,h,H,alpha,f,actual_sol)
    y = [];
    y(1) = alpha;
    w = [];
    w(1) = 0.25;
    w(2) = 0.75;
    k = [];
    t = [a:h:b];
    n = length(t);
    error = [];
    e = 1;
    error(1) = abs(y(1)-actual_sol(a));
    for i=1:n-1
        k(1) = h*f(t(i),y(i));
        k(2) = h*f(t(i) + 2*h/3,y(i) + 2*k(1)/3);
        y(i+1) = y(i) + dot(w,k);
        if mod(i-1,H/h) == 0
            error(e+1) = abs(y(i+1) - actual_sol(t(i+1)));
            e = e+1;
        end
    end
    
    %calculating solution for h=h/2
    h = h/2;
    t = [a:h:b];
    n = length(t);
    y_half = [];
    y_half(1) = alpha;
    e = 1;
    error_half(1) = abs(y(1)-actual_sol(a));
    for i=1:n-1
        k(1) = h*f(t(i),y_half(i));
        k(2) = h*f(t(i) + 2*h/3,y_half(i) + 2*k(1)/3);
        y_half(i+1) = y_half(i) + dot(w,k);
        if mod(i-1,H/h) == 0
        error_half(e+1) = abs(y_half(i+1) - actual_sol(t(i+1)));
        e = e+1;
        end
    end
    
    %Error
     
    fprintf('t\t\tError\t\tError h halved\tratio\n');
    T = [a:H:b];
    for i=1:length(error)
        fprintf('%f\t%f\t%f\t%f\n',T(i),error(i),error_half(i),error_half(i)/error(i));
    end
    
    
end

