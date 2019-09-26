syms x;
f = 2*x/(x^2 - 4);
f_diff = matlabFunction(diff(f));
fprintf('Error bound in (a) = %f\tActual Error = %f\n',Error_Bound(f_diff,1,1.6),Actual_Error(f,1,1.6));

f = exp(3*x)*sin(2*x);
f_diff = matlabFunction(diff(f));
fprintf('Error bound in (b) = %f\tActual Error = %f\n',Error_Bound(f_diff,0,pi/4),Actual_Error(f,0,pi/4));

f = (sin(x))^2 - 2*x*(sin(x)) + 1;
f_diff = matlabFunction(diff(f));
fprintf('Error bound in (c) = %f\tActual Error = %f\n',Error_Bound(f_diff,0.75,1.3),Actual_Error(f,0.75,1.3));

f = 1/(x*log(x));;
f_diff = matlabFunction(diff(f));
fprintf('Error bound in (c) = %f\tActual Error = %f\n',Error_Bound(f_diff,exp(1),exp(1)+1),Actual_Error(f,exp(1),exp(1)+1));

function y = Error_Bound(f,a,b)
    y = 0.5*max(abs(f([a:0.01:b])))*(b-a)^2;
end
function Area = Rectangle(f,a,b)
    Area = f(a)*(b-a);
end
function y = Actual_Error(f,a,b)
g = matlabFunction(f);
Q = integral(g,a,b,'ArrayValued',true);
Q_ = Rectangle(g,a,b);
y = abs(Q-Q_);
end
