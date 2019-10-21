clear all;
clc;
clf;

f = @(t,y) -y + (t.^(0.1)).*(1.1 + t);
actual_sol = @(t) t.^(1.1);
a = 0;
b = 5;
alpha = 0;
h = 0.005;
H = 1;

RK_order2(a,b,h,H,alpha,f,actual_sol);
function RK_order2(a,b,h,H,alpha,f,actual_sol)
    y = [];
    y(1) = alpha;
    y_half = [];
    y_half = alpha;
    h_half = h/2;
    k_half = [];
    w = [];
    w(1) = 0.25;
    w(2) = 0.75;
    k = [];
    t = [a:h:b];
    n = length(t);
    for i=1:n-1
        k(1) = h*f(t(i),y(i));
        k(2) = h*f(t(i) + 2*h/3,y(i) + 2*k(1)/3);
        
        k_half(1) = h_half*f(t(i),y_half(i));
        k_half(2) = h_half*f(t(i) + 2*h_half/3,y_half(i) + 2*k_half(1)/3);
        
        y(i+1) = y(i) + w(1)*k(1) + w(2)*k(2);
        y_half(i+1) = y_half(i) + w(1)*k_half(1) + w(2)*k_half(2);
    end
    T = [a:H:b];
    N = length(T);
    fprintf('t\tError\tError h halved\n');
    for i=1:N
    error = abs(actual_sol(T(i)) - y(int(i*n/N)));
    error_half = abs(actual_sol(T(i)) - y_half(i*n/N));
    fprintf('%f\t%f\t%f\n',t(i),error,error_half);
    end
end