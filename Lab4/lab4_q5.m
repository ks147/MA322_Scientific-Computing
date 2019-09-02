clear all;
clc;
%(a)
x = [0.1 0.2 0.3 0.4];
y = [-0.2900498 -0.56079734 -0.81401972 -1.0526302];

f = matlabFunction(Lagrange(x,y,length(x)));
fprintf('(a) f(0.18) = %e\n\n',f(0.18));

%(b)
x = [-1 -0.5 0 0.5];
y = [0.86199480 0.95802009 1.0986123 1.2943767];
f = matlabFunction(Lagrange(x,y,length(x)));
fprintf('(b) f(0.25 = %e\n\n',f(0.25));
function p = Lagrange(x,y,n)
    syms z;
    p = 0;
    for i=1:n
        prod = y(i);
        for j=1:n
            if j~=i
            prod = prod*((z-x(j))/(x(i)-x(j)));
            end
        end
        p = p+prod;
    end
    
end