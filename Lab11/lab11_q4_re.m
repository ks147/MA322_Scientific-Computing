clear all;
clc;
a = 2;
ic = @(x) 1 + sin(2*pi*x);
bc = @(t) 1;
h = 0.05;
k = 0.01;
a_x = 0;
b_x = 1;
a_t = 0;
b_t = 1;
% [num_sol x t] = Lax_Wendroff(a_x,b_x,a_t,b_t,h,k,a,ic,bc);
% [X,T] = meshgrid(x,t);
% surf(X,T,num_sol');

x = [a_x:h:b_x];
t = [a_t:k:b_t];
m = length(x);
n = length(t);
U = zeros(m,n);
U(:,1) = ic(x);
U(1,:) = bc(t);
U(m,:) = bc(t);
R = a*k/h;
A = zeros(m,m);
%BC
% A(1,1) = 1;
% A(1,2) = -1;
A(m,m) = 1;
for i=2:m-1
    A(i,i-1) = 0.5*(R + R^2);
    A(i,i) = 1- R^2;
    A(i,i+1) = A(i,i-1) - R;
end
for j = 2:n
    U(:,j) = A*U(:,j-1);
end

[X,T] = meshgrid(x,t);
surf(X,T,U');


