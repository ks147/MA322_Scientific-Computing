clear all;
clc;
clf;
close all;
syms x y z 
f1 =  6*x-2*cos(y*z) -1 ;
f2 =  9*y + sqrt(x^2 + sin(z) +1.06) +0.9 ;
f3 =  60*z + 3*exp(-x*y) +10*pi -3 ;


F = jacobian([f1,f2,f3],[x,y,z]);

f1 = matlabFunction(f1);
f2 = matlabFunction(f2);
f3 = matlabFunction(f3);



F11 = matlabFunction(F(1,1));
F12 = matlabFunction(F(1,2));
F13 = matlabFunction(F(1,3));

F21 = matlabFunction(F(2,1));
F22 = matlabFunction(F(2,2));
F23 = matlabFunction(F(2,3));

F31 = matlabFunction(F(3,1));
F32 = matlabFunction(F(3,2));
F33 = matlabFunction(F(3,3));
x= [0;0;0];

xold=x;
f = [f1(0,0,0); f2(0,0,0);f3(0,0,0)];
   
F = [F11(),F12(x(2),x(3)),F13(x(2),x(3));F21(x(1),x(3)),F22(),F23(x(1),x(3));F31(x(1),x(2)),F32(x(1),x(2)),F33()];
x = xold - inv(F)*f;



fprintf('%s %20s %20s %20s\n','n','x','y','z');    
i = 1;
while(norm(x-xold,Inf)>10^(-6))
    
    xold=x;
    f = [f1(x(1),x(2),x(3)); f2(x(1),x(2),x(3));f3(x(1),x(2),x(3))];
    F = [F11(),F12(x(2),x(3)),F13(x(2),x(3));F21(x(1),x(3)),F22(),F23(x(1),x(3));F31(x(1),x(2)),F32(x(1),x(2)),F33()];
    
    fprintf('%d %20f %20f %20f\n',i,x(1),x(2),x(3));    
    x = x - inv(F)*f;
    i = i+1;    
end  