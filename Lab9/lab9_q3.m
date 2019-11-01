clear all;
clc;
%-y'' + p(x)y' + q(x)y = r(x)

p = @(x) 0;
q = @(x) x;
r = @(x) -1;
a = 0;
b = 1;
h = 1/400;
alpha = 1;
beta = 1;
fprintf('Plot for Q3(a)\n\n');
forward_diffQ3a(a,b,h,alpha,beta,p,q,r);

%Q3(b)
p = @(x) 3;
q = @(x) -2;
r = @(x) -2;
a = 0;
b = 1;
h = 1/400;
alpha = 1;
beta = 1;
fprintf('\n\nPlot for Q3(b)\n');
forward_diffQ3b(a,b,h,alpha,beta,p,q,r);


function y = forward_diffQ3a(a,b,h,alpha,beta,p,q,r)
t = [a:h:b];
n = length(t);
y = zeros(1,n);
%Hard coded boundary condition for Q1
A(1,1) = 1-1/h;
A(1,2) = 1/h;
b(1) = beta;
A(n,n) = 1;
b(n,1) = alpha;

for i=2:n-1
    A(i,i-1) =  -1/h^2 - p(t(i))/(2*h);
    A(i,i) = (2/h^2 + q(t(i)));
    A(i,i+1) = (-1/h + p(t(i))/2)/h;
    b(i) = r(t(i));
end

y = A\b;
figure;
plot(t,y);
xlabel('x');
ylabel('y(x)');
title('Second Order Solution Q3(a)');

end
function y = forward_diffQ3b(a,b,h,alpha,beta,p,q,r)
t = [a:h:b];
n = length(t);
y = zeros(1,n);
%Hard coded boundary condition for Q3(b)
A(1,1) = -1-1/h;
A(1,2) = 1/h;
b(1) = alpha;
A(n,n) = 1/h +1;
A(n,n-1) = -1/h;
b(n,1) = beta;

for i=2:n-1
    A(i,i-1) =  -1/h^2 - p(t(i))/(2*h);
    A(i,i) = (2/h^2 + q(t(i)));
    A(i,i+1) = (-1/h + p(t(i))/2)/h;
    b(i) = r(t(i));
end

y = A\b;
figure;
plot(t,y);
xlabel('x');
ylabel('y(x)');
title('Second Order Solution Q3(b)');

end

