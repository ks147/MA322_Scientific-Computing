f = @(x) sin(x);
x = [];
x(1) = 0.3;
h = 0.02;
n = 1+(0.44-0.3)/h;
for i=2:n
    x(i) = x(i-1)+h;
end
y = f(x);
[a b c d] = genspline(x,y,n);
fprintf('The value of spline at 0.3102 = %f\n',evalspline(a,b,c,d,0.3102,x,n));
fprintf('Error at 0.3102 = %e\n',abs(f(0.3102)-evalspline(a,b,c,d,0.3102,x,n)));
function [a b c d] = genspline(x,y,n)
%generates coefficients of cubic spline

for i=1:n-1
    h(i) = x(i+1) - x(i);
end
A = zeros(n,n);
A(1,1) = 1;
A(n,n) = 1;
f = zeros(n,1);

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



