x = [0 3 5 8 13];
y = [0 225 383 623 993];
dy = [75 77 80 74 72];
z = [];
Y = [];

for i=1:length(x)
    z(2*i-1) = x(i);
    z(2*i) = x(i);
    Y(2*i -1) = y(i);
    Y(2*i) = y(i);
end

dist = matlabFunction(Divided_Diff_Interpolation(z,Y,dy,length(z)));
speed = matlabFunction(diff(Divided_Diff_Interpolation(z,Y,dy,length(z))));
acc = matlabFunction(diff(diff(Divided_Diff_Interpolation(z,Y,dy,length(z)))));
fprintf('Position at t=10:%4f,Speed at t=10:%4f\n',dist(10),speed(10));
fplot(speed,[0 12]);
%3(b)
f_55 = @(x) speed(x) - 80.6667;
fprintf('The time at which driver exceeds 55mi/h = %fs\n',Newton(f_55,0.8));
fprintf('Predicted Max speed = %f ft/s\n',speed(Newton(acc,12)));


function p = Divided_Diff_Interpolation(x,y,dy,n)
    
    diff = zeros(n,n);
        
    for i=1:n
        diff(i,1) = y(i);
    end;

    for i= 1:n-1
        for j=0:n-i-1
            if x(j+1)~=x(i+j+1)
                diff(j+1,i+1) = (diff(j+1,i) - diff(j+2,i))/(x(j+1) - x(i+j+1));
            else
                diff(j+1,i+1) = dy(ceil((j+1)/2));
            end
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
function y = Diff(f,x)
    y = (f(x+1e-4)-f(x-1e-4))/(2e-4);
end
    
function root = Newton(f,x0)
x = [];
x(1) = x0;
i = 1;
while(abs(f(x(i))) > 1e-3)
    x(i+1) = x(i) - f(x(i))/Diff(f,x(i));
    i = i+1;
end
root = x(i);
end