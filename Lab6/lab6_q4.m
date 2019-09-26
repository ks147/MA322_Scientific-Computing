%Q4 Compute /integrand 1/(1+4*X^2)
f = @(x) 4/(1+x^2);
a = 0;
b = 1;
actual_integral = integral(f,a,b,'ArrayValued',true);

fprintf('Estimate of pi using\n');
fprintf('(a) Trapezoid Rule = %f\t Error = %f \n',Trapezoid(f,a,b),abs(Trapezoid(f,a,b)-actual_integral));
fprintf('(b) Simpson one-third rule = %f\t Error = %f\n',Simpson(f,a,b),abs(Simpson(f,a,b)-actual_integral));
fprintf('(c) Midpoint rule = %f\t Error = %f\n',Midpt(f,a,b),abs(Midpt(f,a,b)-actual_integral));
fprintf('(d) Simpson three-eigth rule = %f\t Error = %f\n',three_eigth_Simpson(f,a,b),abs(three_eigth_Simpson(f,a,b)-actual_integral));
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
