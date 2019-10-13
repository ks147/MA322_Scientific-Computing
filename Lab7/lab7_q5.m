clear all;
clc;
clf;
close all;
%(a)
f = @(x) (x.^2)*log(x);
a = 1;
b = 1.5;
fprintf('Q5(a)\n');
fprintf('Integral with n = 2: %f\n',gausslege(f,a,b,2));
fprintf('Integral with n = 3: %f\n',gausslege(f,a,b,3));
fprintf('Integral with n = 4: %f\n',gausslege(f,a,b,4));
fprintf('Integral with n = 5: %f\n\n',gausslege(f,a,b,5));

%(b)

f = @(x) exp(3*x)*sin(2*x);
a = 0;
b = pi/4;
fprintf('Q5(b)\n');
fprintf('Integral with n = 2: %f\n',gausslege(f,a,b,2));
fprintf('Integral with n = 3: %f\n',gausslege(f,a,b,3));
fprintf('Integral with n = 4: %f\n',gausslege(f,a,b,4));
fprintf('Integral with n = 5: %f\n\n',gausslege(f,a,b,5));

%(c)
f = @(x) 2/(x.^2 - 4);
a = 0;
b = 0.35;
fprintf('Q5(c)\n');
fprintf('Integral with n = 2: %f\n',gausslege(f,a,b,2));
fprintf('Integral with n = 3: %f\n',gausslege(f,a,b,3));
fprintf('Integral with n = 4: %f\n',gausslege(f,a,b,4));
fprintf('Integral with n = 5: %f\n\n',gausslege(f,a,b,5));

%(d)
f = @(x)  (2.*x)/(x.^2 - 4);
a = 1;
b = 1.6;
fprintf('Q5(d)\n');
fprintf('Integral with n = 2: %f\n',gausslege(f,a,b,2));
fprintf('Integral with n = 3: %f\n',gausslege(f,a,b,3));
fprintf('Integral with n = 4: %f\n',gausslege(f,a,b,4));
fprintf('Integral with n = 5: %f\n\n',gausslege(f,a,b,5));


function I=gausslege(f,a,b,n)
%I=gaussleg(f,a,b,n)
%Aproximates integral using Gauss-Legendre quadrature method
%Legendre polynomial
p=polegende(n);
%Polynomial roots
x=roots(p(n+1,:));
%Change of integration variable if it's needed
if a~=-1 | b~=1
   y=flege(f,a,b);
   G = y(x);
else
   G=f(x);		
end
%derivate
pn=polyder(p(n+1,:));

%Calculus of the coeficients

for i=1:n
   C(i)=2./((1-x(i).^2).*((polyval(pn,x(i))).^2));
    
end

%scalar product of the function at the nodes and the coeficients
I=dot(C,G);
end


function Y=flege(f,a,b)
%Performs variable change if a=!-1 y b=!1
syms x;
x=((b-a)./2).*x+(b+a)./2;
dx=(b-a)./2;
y=feval(f,x)*dx;
Y = matlabFunction(y);
end

function p=polegende(n)
% p=polegend(n)
% Saves on the rows of the p matrix the coeficients of the legendre polin.
p(1,1)=1;
p(2,1:2)=[1 0]; 
for k=2:n
   p(k+1,1:k+1)=((2*k-1)*[p(k,1:k) 0]-(k-1)*[0 0 p(k-1,1:k-1)])/k;
end

end
