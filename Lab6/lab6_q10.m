%length under ellipse
syms x
f = 2*sqrt(1-(x^2)/9);
f_ = matlabFunction(diff(f));
g = @(x) sqrt(1+(f_(x))^2);

[graph_length,n,h] = simpson_error(g,-2.99,2.99,10^-6);
fprintf('length of the graph of the ellipse = %f\n',2*graph_length);

%%functions %%
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

function integral = composite_simpson(f,a,b,n)
h = (b-a)/n;
integral = 0;

for i=1:n
    x = a + (i-1)*h;
    y = a + i*h;
    integral = integral + h*(f(x)+4*f((x+y)/2) + f(y))/6;
end
end

