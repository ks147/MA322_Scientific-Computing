f = @(x) 2*x/(x^2-4);
fprintf('(a) Area under f(x) = %f\n',Rectangle(f,1,1.6));
f = @(x) exp(3*x)*sin(2*x);
fprintf('(b) Area under f(x) = %f\n',Rectangle(f,0,pi/4));
f = @(x) (sin(x))^2 - 2*x*(sin(x)) + 1;
fprintf('(c) Area under f(x) = %f\n',Rectangle(f,0.75,1.3));
f = @(x) 1/(x.*log(x));
fprintf('(d) Area under f(x) = %f\n',Rectangle(f,exp(1),exp(1)+1));

function Area = Rectangle(f,a,b)
    Area = f(a)*(b-a);
end

