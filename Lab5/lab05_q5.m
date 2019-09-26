x = [0.1 0.2 0.3 0.4];
y = [-0.62049958 -0.28398668 0.00660095 0.24842440];
dy = [3.58502082 3.14033271 2.66668043 2.16529366];
n = length(x);
f = @(x) x*cos(x) - 2*x^2 + 3*x -1;

[a b c d] = gen_clamped_spline(x,y,dy,n);

fprintf('The value of spline at 0.2013 = %f\n',evalspline(a,b,c,d,0.2013,x,n));
fprintf('Error at 0.2013 = %e\n',abs(f(0.2013)-evalspline(a,b,c,d,0.2013,x,n)));

function [a b c d] = gen_clamped_spline(x,y,dy,n)
%generates coefficients of cubic spline

for i=1:n-1
    h(i) = x(i+1) - x(i);
end
A = zeros(n,n);
A(1,1) = 2*h(1);
A(1,2) = h(1);
A(n,n) = 2*h(n-1);
A(n,n) = h(n-1);
f = zeros(n,1);
f(1) = 3*((y(2)-y(1))/h(1) - dy(1));
f(n) = 3*(-(y(n)-y(n-1))/h(n-1) + dy(2));
for i=2:n-1
    A(i,i) = 2*(h(i)+h(i-1));
    f(i) = 6*((y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1));
end

for i=2:n-2
    A(i,i+1) = h(i+1);
end

for i=3:n-1
    A(i,i-1) = h(i);
end

%solve for vector s

s = A\f;
%compute coefficients of the spline 
for i=1:n-1
    a(i) = (s(i+1)-s(i))/(6*h(i));
    b(i) = s(i)/2;
    c(i) = (y(i+1)-y(i))/h(i) - (2*h(i)*s(i)+h(i)*s(i+1))/6;
    d(i) = y(i);
end
end

function yy = evalspline(a,b,c,d,xx,x,n)
%evaluate the value of the spline at a point xx
%determine which interval the spline belongs to
i = 1;
while( xx > x(i+1) & i<=n-1)
    i=i+1;
end

%xx is between x(i) and x(i+1)

yy = a(i)*(xx-x(i))^3 + b(i)*(xx-x(i))^2 + c(i)*(xx-x(i)) + d(i);
end
