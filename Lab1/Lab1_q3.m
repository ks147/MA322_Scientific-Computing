a=1;
b=25;
f = @(x) x^3-25

fprintf('Q3 \n');
fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');
i = 0;
while(b-a > 10^-4)
    mid = (a+b)/2;
    fprintf('%2d %15f %15f \n',i,mid,f(mid))
    if f(mid)*f(a) < 0
        b = mid;
    end
    if f(mid)*f(b) < 0
        a = mid;
    end
    if f(mid)==0
        root5 = mid;
        break
    end
    i = i+1;
end
root5 = (a+b)/2;
fprintf('%2d %15f %15f \n',i,root5,f(root5))

