%(a)
f = @(x) 4/(1+x^2);
a = 0;
b = 1;
[ans,sub_interval,actual_error] = error(f,a,b,5*10^-6);
fprintf('Number of sub-intervals to get error < 10^-5/2 = %d\n',ans);
fprintf('Approximate Integral = %f\n\n',composite_simpson(f,a,b,ans));

figure(1);
plot(sub_interval,log(actual_error));
xlabel('Number of Sub intervals');
ylabel('log(Error)');
%(b)
f = @(x) sqrt(1-x^2) - x;
a = 0;
b = 1/sqrt(2);
[ans,sub_interval,actual_error] = error(f,a,b,5*10^-6);
fprintf('Number of sub-intervals to get error < 10^-5/2 = %d\n',ans);
fprintf('Approximate Integral = %f\n\n',composite_simpson(f,a,b,ans));

figure(2);
plot(sub_interval,log(actual_error));
xlabel('Number of Sub intervals');
ylabel('log(Error)');


function [y,sub_interval,actual_error] = error(f,a,b,e)
i = 1;
actual_integ = integral(f,a,b,'ArrayValued',true);
actual_error(1) = abs(composite_simpson(f,a,b,i)-actual_integ);

while(actual_error > e)
    i = i+1;
    actual_error(i) = abs(composite_simpson(f,a,b,i)-actual_integ);
end
y = i;
for i=1:15
    actual_error(i) = abs(composite_simpson(f,a,b,i)-actual_integ);
end
sub_interval = [1:1:15];
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
