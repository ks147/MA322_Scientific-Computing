for i=0:10
    fprintf('The root with constant = %f:%f \n',1.5-0.001*(10-i),Newton(-0.8,10^-3,-1.5-0.001*(10-i)));
end
for i=1:10
    fprintf('The root with constant = %f:%f \n',1.5+0.001*i,Newton(-1,10^-3,-1.5+0.001*i));
end
function y = f_(x,c)
    f = @(x) exp(x)-atan(x)+c;
    h = 10^-5;
    y = (f(x+h)-f(x-h))/(2*h);
end
function y=Newton(a,limit,c)
    f = @(x) exp(x)-atan(x)+c;
               
    x=[];
    x(1) = a;
    i=1;
    while(abs(f(x(i))) > limit)
        x(i+1) = x(i) - (f(x(i))/f_(x(i),c));
        i=i+1;
    end
    y = x(i);
end
