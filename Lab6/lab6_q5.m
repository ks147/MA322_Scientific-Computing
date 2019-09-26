
%(a)
f = @(x) x^2*exp(-x^2);
fprintf('Estimate of integral = %f\tActual Error = %f\n',composite_trapezoid(f,0,2,8),Actual_Error(f,0,2,8));
%(b)
f = @(x) 1/(x*log(x));
fprintf('Estimate of integral = %f\tActual Error = %f\n',composite_trapezoid(f,exp(1),exp(1)+2,8),Actual_Error(f,exp(1),exp(1)+2,8));
%(c)
f = @(x) x^2*log(x^2+1);
fprintf('Estimate of integral = %f\tActual Error = %f\n',composite_trapezoid(f,0,2,8),Actual_Error(f,0,2,8));
%(d)
f = @(x) (sin(x))^2 - 2*x*sin(x) + 1;
fprintf('Estimate of integral = %f\tActual Error = %f\n',composite_trapezoid(f,0.75,1.75,8),Actual_Error(f,0.75,1.75,8));


function integral = composite_trapezoid(f,a,b,n)
h = (b-a)/n;
integral = 0;
for i=1:n
    integral = integral + h*(f(a+(i-1)*h) + f(a+i*h))/2;
end
end

function y = Actual_Error(f,a,b,n)

Q = integral(f,a,b,'ArrayValued',true);
Q_ = composite_trapezoid(f,a,b,n);
y = abs(Q-Q_);
end