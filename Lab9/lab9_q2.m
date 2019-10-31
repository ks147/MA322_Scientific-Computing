
function y = forward_diff(a,b,h,alpha,beta,p,q,r,fig_no)
t = [a:h:b];
n = length(t);
y = zeros(1,n);
y(1) = alpha;
y(n) = beta;

A = zeros(n-2,n-2);
b = zeros(n-2,1);
% i = 2
A(1,1) = (2/h^2 - p(t(2))/h + q(t(2)));
A(1,2) = (-1/h + p(t(2)))/h;
b(1) = r(t(2)) - y(1)*(-1/h^2);
%i = n
A(n-2,n-3) = -1/h^2;
A(n-2,n-2) = (2/h^2 - p(t(n-1))/h + q(t(n-1)));
b(n-2) = r(t(n-1)) - y(n)*(-1/h + p(t(n-1)))/h;

for i=3:n-2
    A(i-1,i-2) =  -1/h^2;
    A(i-1,i-1) = (2/h^2 - p(t(i))/h + q(t(i)));
    A(i-1,i) = (-1/h + p(t(i)))/h;
    b(i-1) = r(t(i));
end
y(2:n-1) = A\b;
figure(fig_no);
plot(t,y);
xlabel('x');
ylabel('y(x)');
title('Forward Difference');

end



