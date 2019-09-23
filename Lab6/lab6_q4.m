%Q4 Compute /integrand 1/(1+4*X^2)
f = @(x) 4/(1+x^2);
a = 0;
b = 1;

fprintf('Estimate of pi using\n');
fprintf('(a) Trapezoid Rule = %f\n',Trapezoid(f,a,b));
fprintf('(b) Simpson one-third rule = %f\n',Simpson(f,a,b));
fprintf('(c) Midpoint rule = %f\n',Midpt(f,a,b));
fprintf('(d) Simpson three-eigth rule = %f\n',three_eigth_Simpson(f,a,b));
function y = Trapezoid(f,a,b)
    y = 0.5*(f(a)+f(b))*(b-a);
end
function y = Simpson(f,a,b)
    y = (b-a)*(f(a)+4*f((a+b)/2)+f(a))/6;
end
function y = Midpt(f,a,b)
    y = (b-a)*(f((a+b)/2));
end
function y = three_eigth_Simpson(f,a,b)
    h = (b-a)/3;
    y = 3*h*(f(a)+3*f(a+h)+3*f(a+2*h)+f(b))/8;
end
