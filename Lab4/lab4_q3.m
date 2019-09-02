clear all;
clc;

f = @(x) log(x);
x = [1 1.1 1.3 1.4];
Lagrange(f,x,4);

function Lagrange(f,x,n)
    syms z;
    p = 0;
    for i=1:n
        prod = f(x(i));
        for j=1:n
            if j~=i
            prod = prod*((z-x(j))/(x(i)-x(j)));
            end
        end
        p = p+prod;
    end
    F(z) = f(z);
    Error = matlabFunction(F - p);
    E = @(x) abs(Error(x));
    fplot(E,[1 1.4]);
end