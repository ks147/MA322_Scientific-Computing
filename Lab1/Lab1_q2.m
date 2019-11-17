%q2(a)
a=0.5;
b=1.5;

fprintf('Q2(a) \n');
fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');

f = @(x) 2+cos(exp(x)-2)-exp(x)
i = 0;
while(b-a > 10^-3)
    mid = (a+b)/2;
    fprintf('%2d %15f %15f \n',i,mid,f(mid))
    if f(mid)*f(a) < 0
        b = mid;
    end
    if f(mid)*f(b) < 0
        a = mid;
    end
    if f(mid)==0
        root1 = mid;
        break
    end
    i=i+1;
end

root1 = (a+b)/2;

%q2(b)
a=0;
b=1;
f = @(x) x-2^(-x)
fprintf('Q2(b) \n');
fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');
while(b-a > 10^-3)
    mid = (a+b)/2;
    fprintf('%2d %15f %15f \n',i,mid,f(mid))
    if f(mid)*f(a) < 0
        b = mid;
    end
    if f(mid)*f(b) < 0
        a = mid;
    end
    if f(mid)==0
        root2 = mid;
        break
    end
    i=i+1;
end

root2 = (a+b)/2;

%q2(c)
a=0;
b=1;

f = @(x) exp(x)-x^2+3*x-2
fprintf('Q2(c) \n');
fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');
i=0;
while(b-a > 10^-3)
    mid = (a+b)/2;
    fprintf('%2d %15f %15f \n',i,mid,f(mid))
    if f(mid)*f(a) < 0
        b = mid;
    end
    if f(mid)*f(b) < 0
        a = mid;
    end
    if f(mid)==0
        root3 = mid;
        break
    end
    i=i+1;
end

root3 = (a+b)/2;

%q2(d)
a=-3;
b=-2;
fprintf('Q2(d) \n');
fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');
f = @(x) 2*x*cos(2*x)-(x+1)^2
i=0;
while(b-a > 10^-3)
    mid = (a+b)/2;
    fprintf('%2d %15f %15f \n',i,mid,f(mid))
    if f(mid)*f(a) < 0
        b = mid;
    end
    if f(mid)*f(b) < 0
        a = mid;
    end
    if f(mid)==0
        root3 = mid;
        break
    end
    i=i+1;
end

root3 = (a+b)/2;

a=-1;
b=0;
i=0;
fprintf('Q2(d) \n');
fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');

while(b-a > 10^-3)
    mid = (a+b)/2;
    fprintf('%2d %15f %15f \n',i,mid,f(mid))
    if f(mid)*f(a) < 0
        b = mid;
    end
    if f(mid)*f(b) < 0
        a = mid;
    end
    if f(mid)==0
        root4 = mid;
        break
    end
    i=i+1;
end


root4 = (a+b)/2;

%q2(e)
a=1.2;
b=1.3;
f = @(x) x*cos(x)-2*x^2+3*x-1
fprintf('Q2(e) \n');
fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');
i = 0;
while(b-a > 10^-3)
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
    i=i+1;
end
root5 = (a+b)/2;


