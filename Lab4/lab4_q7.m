clear all;
clc;
x = [2,3,5,6];
y = [1.5713,1.5719,1.5738,1.5751];
f = matlabFunction(Lagrange(x,y,3));
fprintf('Interpolating using second degree polynomial f(4) = %.4e\n\n',f(4));

f = matlabFunction(Lagrange(x,y,4));
fprintf('Interpolating using third degree polynomial f(4) = %.4e\n\n',f(4));

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