clear all;
clc;
x = [0 0.1 0.3 0.6 1.0];
y = [-6 -5.89483 -5.65014 -5.17788 -4.28172];
p1 = matlabFunction(Divided_Diff_Interpolation(x,y,5));
p2 = matlabFunction(Lagrange(x,y,5));
fprintf('Value of f(0.2) using Divided Diff Interpolation:%e\n',p1(0.2));
fprintf('Value of f(0.2) using Lagrange Interpolation:%e\n\n',p2(0.2));
x(6) = 1.1;
y(6) = -3.99583;
p1 = matlabFunction(Divided_Diff_Interpolation(x,y,6));
p2 = matlabFunction(Lagrange(x,y,6));
fprintf('New Value of f(0.2) using Divided Diff Interpolation:%e\n',p1(0.2));
fprintf('New Value of f(0.2) using Lagrange Interpolation:%e\n',p2(0.2));
function p = Divided_Diff_Interpolation(x,y,n)
    
    diff = zeros(n,n);
        
    for i=1:n
        diff(i,1) = y(i);
    end;

    for i= 1:n-1
        for j=0:n-i-1
            diff(j+1,i+1) = (diff(j+1,i) - diff(j+2,i))/(x(j+1) - x(i+j+1));
        end
    end
%y[j][i] = (y[j][i - 1] - y[j + 1][i - 1]) / (x[j] - x[i + j]); 
    syms z
   
    p = y(1,1);
    
    for i = 2:n
        p = p + (u_cal_for(z,x,i)*diff(1,i));
    end;
   
end

function y = u_cal_for(u,x,n)
    y = 1;
    for i=1:n-1
        y = y*(u-x(i));
    end
end
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