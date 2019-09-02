
x = [1950:10:2000];
y = [151326 179323 203302 226542 249633 281422];
f = matlabFunction(Divided_Diff_Interpolation(x,y,length(x)));
fprintf('Population in 1940 = %d\n',f(1940));
fprintf('Population in 1975 = %d\n',f(1975));
fprintf('Population in 2020 = %d\n',f(2020));

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