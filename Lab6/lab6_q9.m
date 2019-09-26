f = @(x) x*log(x);
a = 1;
b = 2;
e = 10^-5;
%Composite Trapezoid
[appx_integral,n,h] = trapezoid_error(f,a,b,e);
fprintf('Approximation using Composite Trapezoid = %f\tn = %f\th = %f\n',appx_integral,n,h);

%Composite Simpson
[appx_integral,n,h] = simpson_error(f,a,b,e);
fprintf('Approximation using Composite Simpson = %f\tn = %f\th = %f\n',appx_integral,n,h);

%Composite Midpoint 
[appx_integral,n,h] = midpt_error(f,a,b,e);
fprintf('Approximation using Composite Midpoint = %f\tn = %f\th = %f\n',appx_integral,n,h);

%% functions %%

function [appx_integral,n,h] = simpson_error(f,a,b,e)
i = 1;
actual_integ = integral(f,a,b,'ArrayValued',true);
actual_error(1) = abs(composite_simpson(f,a,b,i)-actual_integ);

while(actual_error > e)
    i = i+1;
    actual_error(i) = abs(composite_simpson(f,a,b,i)-actual_integ);
end
n = i;
h = (b-a)/n;
appx_integral = composite_simpson(f,a,b,n);

end

function [appx_integral,n,h] = midpt_error(f,a,b,e)
i = 1;
actual_integ = integral(f,a,b,'ArrayValued',true);
actual_error(1) = abs(composite_midpt(f,a,b,i)-actual_integ);

while(actual_error > e)
    i = i+1;
    actual_error(i) = abs(composite_midpt(f,a,b,i)-actual_integ);
end
n = i;
h = (b-a)/n;
appx_integral = composite_midpt(f,a,b,n);

end

function [appx_integral,n,h] = trapezoid_error(f,a,b,e)
i = 1;
actual_integ = integral(f,a,b,'ArrayValued',true);
actual_error(1) = abs(composite_trapezoid(f,a,b,i)-actual_integ);

while(actual_error > e)
    i = i+1;
    actual_error(i) = abs(composite_trapezoid(f,a,b,i)-actual_integ);
end
n = i;
h = (b-a)/n;
appx_integral = composite_trapezoid(f,a,b,n);

end

function integral = composite_trapezoid(f,a,b,n)
h = (b-a)/n;
integral = 0;
for i=1:n
    integral = integral + h*(f(a+(i-1)*h) + f(a+i*h))/2;
end

end

function integral = composite_midpt(f,a,b,n)
h = (b-a)/n;
integral = 0;
for i=1:n
    x = a + (i-1)*h;
    y = a + i*h;
    integral = integral + h*f((x+y)/2);

end
end

function integral = composite_simpson(f,a,b,n)
h = (b-a)/n;
integral = 0;

for i=1:n
    x = a + (i-1)*h;
    y = a + i*h;
    integral = integral + h*(f(x)+4*f((x+y)/2) + f(y))/6;
end
end

