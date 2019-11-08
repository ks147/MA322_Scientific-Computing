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
h = 0.04;
k = h^2/4;                                 %k = h^2/4
c = 4/(pi^2);
exact_sol = @(x,t) exp(-t).*sin(pi*x/2) + exp(-t/4).*sin(pi*x/4);

% fprintf('\t\t FTCS\n');
% [num_sol x t] = ftcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
% plot_end_time(exact_sol,num_sol,x,t);
% surface_plot(exact_sol,num_sol,x,t);
% % plot_N_Max_error('ftcs',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);
% 
% fprintf('\t\tBTCS\n');
% [num_sol x t] = btcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
% plot_end_time(exact_sol,num_sol,x,t);
% surface_plot(exact_sol,num_sol,x,t);
% % plot_N_Max_error('btcs',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);
% 
% fprintf('\t\tCrank Nicolson\n');
% [num_sol x t] = Crank_Nicolson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
% plot_end_time(exact_sol,num_sol,x,t);
% surface_plot(exact_sol,num_sol,x,t);
% %plot_N_Max_error('crank nicolson',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);
% 
% fprintf('\t\tRichardson\n');
% [num_sol x t] = Richardson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
% plot_end_time(exact_sol,num_sol,x,t);
% surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('richardson',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);

%%(b)
ic = @(x) sin(pi.*x);      %Initial Condition
bc_1 = @(t) 0;                                  %Boundary Condition 1
bc_2 = @(t) 0;                                  %Boundary Condition 2
a_x = 0;                                        %Range x 
b_x = 1;
a_t = 0;                                        %Range t
b_t = 1;
h = 0.04;
k = h^2/4;                                 %k = h^2/4
c = 1;
exact_sol = @(x,t) exp(-(pi*pi*t)).*sin(pi*x);
fprintf('\t\t FTCS\n');
[num_sol x t] = ftcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('ftcs',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);

fprintf('\t\tBTCS\n');
[num_sol x t] = btcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('btcs',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);

fprintf('\t\tCrank Nicolson\n');
[num_sol x t] = Crank_Nicolson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('crank nicolson',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);

fprintf('\t\tRichardson\n');
[num_sol x t] = Richardson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('richardson',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);


%%(c)%%%%
ic = @(x) sin(pi.*x/2) + 0.5*sin(2*pi*x);      %Initial Condition
bc_1 = @(t) 0;                                  %Boundary Condition 1
bc_2 = @(t) exp(-(pi*pi*t/4));                                  %Boundary Condition 2
a_x = 0;                                        %Range x 
b_x = 1;
a_t = 0;                                        %Range t
b_t = 1;
h = 0.04;
k = h^2/4;                                 %k = h^2/4
c = 1;
exact_sol = @(x,t) exp(-(pi*pi*t/4)).*sin(pi*x/2) + 0.5*exp(-(pi*pi*t*4)).*sin(pi*x*2) ;
fprintf('\t\t FTCS\n');
[num_sol x t] = ftcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('ftcs',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);

fprintf('\t\tBTCS\n');
[num_sol x t] = btcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('btcs',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);

fprintf('\t\tCrank Nicolson\n');
[num_sol x t] = Crank_Nicolson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('crank nicolson',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);

fprintf('\t\tRichardson\n');
[num_sol x t] = Richardson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);
% plot_N_Max_error('richardson',c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol);


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
    U(i,m) = bc_2(t(i));
    for j=2:m-1
        U(i,j) = s*U(i-1,j-1) + (1 - 2*s)*U(i-1,j) + s*U(i-1,j+1);
    end
end

end

%implicit method
function [U x t] = btcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2)
s = c*(k/(h^2));            %lambda = k/h^2;

x = [a_x:h:b_x];            %Discretize space
t = [a_t:k:b_t];            %Discretize time range
n = length(t);
m = length(x);
U = zeros(n,m);
U(1,:) = ic(x);
%U(t,x)
%AU(j,:) = AU(j-1,:)
A = zeros(m,m);
A(1,1) = 1+2*s;
A(1,2) = -s;
A(m,m-1) = -s;
A(m,m) = 1+2*s;
for i = 2:m-1
    A(i,i-1) = -s;
    A(i,i) = 1+2*s;
    A(i,i+1) =  -s;
end
U = U';
for j = 2:n
   U(1,j-1) = bc_1(t(j-1));
   U(m,j-1) = bc_2(t(j-1));
   U(:,j) = A\U(:,j-1); 
 
end
U(1,n) = bc_1(t(n));
U(m,n) = bc_2(t(n));
U = U';
end

function [U x t] = Crank_Nicolson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2)
s = c*(k/(h^2));            %lambda = k/h^2;

x = [a_x:h:b_x];            %Discretize space
t = [a_t:k:b_t];            %Discretize time range
n = length(t);
m = length(x);
U = zeros(n,m);
U(1,:) = ic(x);
%U(t,x)
%(2I + sB)U(j,:) = (2I-sB)U(j-1,:)
B = zeros(m,m);
B(1,1) = 2;
B(1,2) = -1;
B(m,m-1) = -1;
B(m,m) = 2;
for i = 2:m-1
    B(i,i-1) = -1;
    B(i,i) = 2;
    B(i,i+1) = -1;
end
U = U';
for j = 2:n
   U(1,j-1) = bc_1(t(j-1));
   U(m,j-1) = bc_2(t(j-1));
   U(:,j) = (2*eye(m,m) + s*B)\(2*eye(m,m) - s*B)*U(:,j-1); 
end
U(1,n) = bc_1(t(n));
U(m,n) = bc_2(t(n));
U = U';
end

%Explicit scheme hence not solved using matrix
%Can construct matrix like in btcs,crank_nicolson
function [U x t] = Richardson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2)
s = c*(k/(h^2));          %lambda = k/h^2;

x = [a_x:h:b_x];    %Discretize space
t = [a_t:k:b_t];    %Discretize time range
n = length(t);
m = length(x);
U = zeros(n,m);
U(1,:) = ic(x);
for i=2:n
    U(i,1) = bc_1(t(i));
    U(i,m) = bc_2(t(i));
    for j=2:m-1
        if i==2
        U(i,j) = s*U(i-1,j-1) + (1 - 2*s)*U(i-1,j) + s*U(i-1,j+1);
        else
        U(i,j) = (1-2*s)*U(i-2,j) + 2*s*U(i-1,j+1) + 2*s*U(i-1,j-1);
        U(i,j) = U(i,j)/(1+2*s);
        end
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
legend('Numerical Solution','Exact Solution');
title('Plot of exact and Numerical Solution at Final time');


end

function surface_plot(exact_sol,num_sol,x,t)
    figure;
    subplot(1,2,1);
    mesh(x,t,num_sol);
    subplot(1,2,2);
    [X Y] = meshgrid(x,t);
    surf(X,Y,exact_sol(X,Y));
    title('Surface plot of Numerical Solution and Exact solution');
end

%scheme-ftcs,btcs,Crank Nicolson etc.
function plot_N_Max_error(scheme,c,a_x,b_x,a_t,b_t,ic,bc_1,bc_2,exact_sol) 
    i = 1;
    N = 1;
    for N = 10:20:90
        h = (b_x - a_x)/(N-1);
        k = h^2/8;
        if scheme=='ftcs'
            [U x t] = ftcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
        elseif scheme=='btcs'
            [U x t] = btcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
        elseif scheme=='crank nicolson'
            [U x t] = Crank_Nicolson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
        else
            [U x t] = Richardson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
        end
        error = abs(exact_sol(x',t) - U');
        Max_error(i) = max(max(error));
        i = i+1;
        
    end
    figure;
    plot([10:20:90],log(Max_error));
    xlabel('N');
    ylabel('log(Max Error)');
    title('N vs Error');    

end
