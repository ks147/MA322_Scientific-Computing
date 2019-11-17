clear all;
clc;

%du/dt = (2/pi)^2 d^2u/dx^2
%y = pi/2x

ic = @(x) sin(pi.*x/4).*(1 + 2.*cos(pi*x/4));      %Initial Condition
bc_1 = @(t) 0;                                  %Boundary Condition 1
bc_2 = @(t) 0;                                  %Boundary Condition 2
a_x = 0;                                        %Range x 
b_x = 4;
a_t = 0;                                        %Range t
b_t = 1;
h = 0.01;
k = h^2/6;                                 %k = h^2/4
c = 4/(pi^2);
exact_sol = @(x,t) exp(-t).*sin(pi*x/2) + exp(-t/4).*sin(pi*x/4);

[num_sol x t] = Duford(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);

%%FUNCTIONS
function [U x t] = ftcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2)
    
s = c*(k/(h^2));          %lambda = k/h^2;

x = [a_x:h:b_x];    %Discretize space
t = [a_t:k:b_t];    %Discretize time range
n = length(t);
m = length(x);
U = zeros(n,m);
U(1,:) = ic(x);
for i=2:n
    U(i,1) = bc_1(t(i));
    U(i,m-1) = bc_2(t(i));
    for j=2:m-1
        U(i,j) = s*U(i-1,j-1) + (1 - 2*s)*U(i-1,j) + s*U(i-1,j+1);
    end
end

end

function plot_end_time(exact_sol,num_sol,x,t)
figure;
[n m] = size(num_sol);
plot(x,num_sol(n,:),'g+');
hold on;
plot(x,exact_sol(x,t(n)),'b');
hold off;
legend('Numerical Solution','Exact Solution')

end
function surface_plot(exact_sol,num_sol,x,t)
figure;
[X,T] = meshgrid(x,t);
surf(X,T,num_sol);
figure;
surf(X,T,exact_sol(X,T));
end



