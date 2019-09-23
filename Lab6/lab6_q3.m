
for i=1:4 
if i==1    
f = @(x) 2*x/(x^2-4);
a = 1;
b = 1.6;
elseif i==2
    f = @(x) exp(3*x)*sin(2*x);
    a = 0;
    b = pi/4;
elseif i==3
    f = @(x) exp(3*x)*sin(2*x);
    a = 0;
    b = pi/4;
else 
    f = @(x) 1/(x.*log(x));
    a = exp(1);
    b = exp(2);
end

    fprintf('(%d) Area under f(x)\n\tEstimate by Midpt = %f\n',i,Midpt(f,a,b));
    fprintf('\tEstimate by Trapezoid Rule = %f\n',Trapezoid(f,a,b));
    fprintf('\tEstimate by Simpson Rule = %f\n',Simpson(f,a,b));
    fprintf('\tEstimate by Corrected Trapezoid Rule = %f\n',Corr_Trapezoid(f,a,b));
    fprintf('\tActual Area = %f\n',integral(f,a,b,'ArrayValued',true));
    
end

function y = Trapezoid(f,a,b)
    y = 0.5*(f(a)+f(b))*(b-a);
end
function y = Simpson(f,a,b)
    y = (b-a)*(f(a)+4*f((a+b)/2)+f(a))/6;
end
function y = Midpt(f,a,b)
    y = (b-a)*(f((a+b)/2));
end
function y = Corr_Trapezoid(f,a,b)
    y = (b-a)*(f(a)+f(b))/2 + (b-a)^2*(f(a)-f(b))/12;
end

