%(a)
f = @(x) x^2*exp(-x^2);
fprintf('Estimate using composite Midpoint = %f\t using composite Simpson = %f\n',composite_midpt(f,0,2,8),composite_simpson(f,0,2,8));
%(b)
f = @(x) 1/(x*log(x));
fprintf('Estimate using composite Midpoint = %f\t using composite Simpson = %f\n',composite_midpt(f,exp(1),exp(1)+2,8),composite_simpson(f,exp(1),exp(1)+2,8));
%(c)
f = @(x) x^2*log(x^2+1);
fprintf('Estimate using composite Midpoint = %f\t using composite Simpson = %f\n',composite_midpt(f,0,2,8),composite_simpson(f,0,2,8));
%(d)
f = @(x) (sin(x))^2 - 2*x*sin(x) + 1;
fprintf('Estimate using composite Midpoint = %f\t using composite Simpson = %f\n',composite_midpt(f,0.75,1.75,8),composite_simpson(f,0.75,1.75,8));
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

