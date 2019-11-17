%8(a)
fprintf('8(a) has infinite roots \n')
%8(b)
fprintf('8(b) has inifinte roots \n');
%8(c)
fprintf('8(c) has roots: \n');
Secant(0.3,0.5,10^-5,'c');
Secant(-1,0,10^-5,'c');
Secant(10,12,10^-5,'c');

%8(d)
fprintf('\n8(d) has roots: \n');
Secant(-0.6,0,10^-5,'d');
Secant(-1,-0.7,10^-5,'d');
Secant(-3,-2,10^-5,'d');

%8(e)
fprintf('\n8(e) has roots: \n');
Secant(-2,-1,10^-5,'e');
Secant(1,2,10^-5,'e');

function Secant(a,b,limit,part)

    if(part=='a')
        f = @(x) sin(x/2) - 1;
    elseif(part=='b')
        f = @(x) exp(x) - tan(x);
    elseif(part=='c')
        f = @(x) x^3 - 12*(x^2) + 3*x+1;
    elseif(part=='d')
        f = @(x) x^3 + 4.001*(x^2) + 4.002*x + 1.101;
    elseif(part=='e')
        f = @(x) x^6 - x^4 + 2*(x^3) - 3*(x^2) + x -4 ;
    end
    x=[];
    x(1) = b;
    x(2) = a;
    i=2;
    fprintf('%f <= x <= %f \n',a,b);
    fprintf('%2s %15s %15s \n','n','x(n)','f(x(n))');

    while(abs(f(x(i))) > limit)
       
        derivative = f(x(i))-f(x(i-1));
        derivative = derivative/(x(i)-x(i-1));
        x(i+1) = x(i) - (f(x(i))/derivative);
        i=i+1;
    end
    fprintf('%2d %15f %15f \n',i,x(i),f(x(i)));
end