%q1
x=zeros(1,1000);
x(1)=-4;
x(2)=-3;
a=x(1);
b=x(2);
ans=-2;

fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');
for i=3:100
    x(i) = (a+b)/2;
    if f(x(i))*f(a) < 0
        b = x(i);
    end
    if f(x(i))*f(b) < 0
        a = x(i);
    end
    if f(x(i))==0
        ans = x(i);
        break
    end
    if(mod(i,5)==0)
       fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
    end
end
if ans==-2
ans = x(1000);
end


function y=f(x)
    y = exp(x)-sin(x);
end


