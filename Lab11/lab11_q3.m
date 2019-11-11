clc;
close all;

%1D Hyperbolic Equation
a = -2;
a_x = 0;
b_x = 1;
a_t = 0;
b_t = 1;
ic = @(x) 1 + sin(2*pi*x);
bc = @(t) 1;
h = 0.005;
k = 0.001;
[num_sol x t] = ftfs_1d_hyperbolic(a,a_x,b_x,a_t,b_t,h,k,ic,bc);
figure;
plot(x,num_sol(:,3));
xlabel('x')
ylabel('U(x,0.003)');
title('Numerical Solution Using FTFS');
[num_sol x t] = btfs_1d_hyperbolic(a,a_x,b_x,a_t,b_t,h,k,ic,bc);
figure;
plot(x,num_sol(:,3));
xlabel('x')
ylabel('U(x,0.003)');
title('Numerical Solution Using BTFS');