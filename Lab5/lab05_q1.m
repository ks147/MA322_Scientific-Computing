
x = [-1 -0.5 0 0.5];
y = [0.86199480 0.95802009 1.0986123 1.2943767];
dy = [0.15536240 0.23269654 0.33333333 0.45186776];
z = [];
Y = [];

for i=1:length(x)
    z(2*i-1) = x(i);
    z(2*i) = x(i);
    Y(2*i -1) = y(i);
    Y(2*i) = y(i);
end
f = @(x) log(exp(x)+2);
p = matlabFunction(Divided_Diff_Interpolation(z,Y,dy,length(z)));
fprintf('Interpolated Value f(0.25) = %e\n',p(0.25));
fprintf('Absolute Error = %e\n',abs(p(0.25)-f(0.25)));
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