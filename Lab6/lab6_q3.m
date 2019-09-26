clear all;
clc;

%Part a
syms x;
df = (2*x)/(x^2 - 4);
df = diff(df);
ddf = diff(df);
df3 = diff(ddf);
df4 = diff(df3);
df = matlabFunction(df);
ddf = matlabFunction(ddf);
df4 = matlabFunction(df4);
f = @(x) (2.*x)./(x.^2 - 4);
a = 1;
b = 1.6;

range = a:0.01:b;
mxddf = max(abs(ddf(range)));

I = integral(f,a,b);
val = (b-a)*f((a+b)/2);
fprintf('By Midpoint rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/24;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a)+f(b));
fprintf('By Trapezoid rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/12;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

mxdf4 = max(abs(df4(range)));
val = ((b-a)/6)*(f(a) + 4*f((a+b)/2) + f(b));
fprintf('By Simpsons rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*((b-a)/2)^5)/90;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a) + f(b)) - ((((b-a)^2)/12)*(df(b) - df(a)));
fprintf('By Corrected Trapezoidal rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*(b-a)^5)/720;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n\n', actual_error);

%Part b
syms x;
df = exp(3*x)*sin(2*x);
df = diff(df);
ddf = diff(df);
df3 = diff(ddf);
df4 = diff(df3);
df = matlabFunction(df);
ddf = matlabFunction(ddf);
df4 = matlabFunction(df4);
f = @(x) exp(3.*x).*sin(2.*x);
a = 0;
b = pi/4;

range = a:0.01:b;
mxddf = max(abs(ddf(range)));

I = integral(f,a,b);
val = (b-a)*f((a+b)/2);
fprintf('By Midpoint rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/24;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a)+f(b));
fprintf('By Trapezoid rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/12;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

mxdf4 = max(abs(df4(range)));
val = ((b-a)/6)*(f(a) + 4*f((a+b)/2) + f(b));
fprintf('By Simpsons rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*((b-a)/2)^5)/90;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a) + f(b)) - ((((b-a)^2)/12)*(df(b) - df(a)));
fprintf('By Corrected Trapezoidal rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*(b-a)^5)/720;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n\n', actual_error);

%Part c
syms x;
df = (2*x)/(x^2 - 4);
df = diff(df);
ddf = diff(df);
df3 = diff(ddf);
df4 = diff(df3);
df = matlabFunction(df);
ddf = matlabFunction(ddf);
df4 = matlabFunction(df4);
f = @(x) (sin(x).*sin(x) - 2.*x.*sin(x) + 1);
a = 0.75;
b = 1.3;

range = a:0.01:b;
mxddf = max(abs(ddf(range)));

I = integral(f,a,b);
val = (b-a)*f((a+b)/2);
fprintf('By Midpoint rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/24;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a)+f(b));
fprintf('By Trapezoid rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/12;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

mxdf4 = max(abs(df4(range)));
val = ((b-a)/6)*(f(a) + 4*f((a+b)/2) + f(b));
fprintf('By Simpsons rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*((b-a)/2)^5)/90;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a) + f(b)) - ((((b-a)^2)/12)*(df(b) - df(a)));
fprintf('By Corrected Trapezoidal rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*(b-a)^5)/720;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n\n', actual_error);

%Part d
syms x;
df = (2*x)/(x^2 - 4);
df = diff(df);
ddf = diff(df);
df3 = diff(ddf);
df4 = diff(df3);
df = matlabFunction(df);
ddf = matlabFunction(ddf);
df4 = matlabFunction(df4);
f = @(x) 1./(x.*log(x));
a = exp(1);
b = exp(1) + 1;

range = a:0.01:b;
mxddf = max(abs(ddf(range)));

I = integral(f,a,b);
val = (b-a)*f((a+b)/2);
fprintf('By Midpoint rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/24;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a)+f(b));
fprintf('By Trapezoid rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxddf*(b-a)^3)/12;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

mxdf4 = max(abs(df4(range)));
val = ((b-a)/6)*(f(a) + 4*f((a+b)/2) + f(b));
fprintf('By Simpsons rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*((b-a)/2)^5)/90;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n', actual_error);

val = ((b-a)/2)*(f(a) + f(b)) - ((((b-a)^2)/12)*(df(b) - df(a)));
fprintf('By Corrected Trapezoidal rule, estimated value is %.5f.\n', val);
actual_error = abs(I-val);
error_bound = (mxdf4*(b-a)^5)/720;
fprintf('Bound for error by error formula is %.5f.\n', error_bound);
fprintf('Actual error is %.5f.\n\n', actual_error);