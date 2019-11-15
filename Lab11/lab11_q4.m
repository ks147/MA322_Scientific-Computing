clear all;
clc;

h=0.05;
k=0.01;
xa=0;
xb=1;
ta=0;
tb=1;
x=[xa:h:xb];
t=[ta:k:tb];
nx=(xb-xa)/h+1; %no of nodal points in x-axis
nt=(tb-ta)/k+1;
lambda=k/h;
a=2;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(x);
u(1,:)=bc_xa(t);
u(nx,:)=bc_xb(t);

for j=2:nt
for i=2:nx-1
        u(i,j)=u(i,j-1) - a*lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*a*a*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2; ;
end
end
[X,Y]=meshgrid(t,x);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('(i)Lax Wendroff Approx.');


%b
clear all;
clc;

h=0.05;
k=0.01;
xa=0;
xb=1;
ta=0;
tb=1;
x=[xa:h:xb];
t=[ta:k:tb];
nx=(xb-xa)/h+1; %no of nodal points in x-axis
nt=(tb-ta)/k+1;
lambda=k/h;
a=2;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(x);
u(1,1)=1;
u(nx,:)=bc_xb(t);

for j=2:nt
    for i=2:nx-1
        u(i,j)=u(i,j-1) - a*lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*a*a*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2; 
    end
    u(1,j)=u(2,j);
end
[X,Y]=meshgrid(t,x);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('(ii)Lax Wendroff Approx.');


%c
clear all;
clc;

h=0.05;
k=0.01;
xa=0;
xb=1;
ta=0;
tb=1;
x=[xa:h:xb];
t=[ta:k:tb];
nx=(xb-xa)/h+1; %no of nodal points in x-axis
nt=(tb-ta)/k+1;
lambda=k/h;
a=2;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(x);
u(1,1)=1;
u(nx,:)=bc_xb(t);

for j=2:nt
    for i=2:nx-1
        u(i,j)=u(i,j-1) - a*lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*a*a*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2;
    end
    u(1,j)=u(1,j-1) - a*lambda*(u(2,j-1) - u(1,j-1));
end
[X,Y]=meshgrid(t,x);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('(iii)Lax Wendroff Approx.');


%d
clear all;
clc;

h=0.05;
k=0.01;
xa=0;
xb=1;
ta=0;
tb=1;
x=[xa:h:xb];
t=[ta:k:tb];
nx=(xb-xa)/h+1; %no of nodal points in x-axis
nt=(tb-ta)/k+1;
lambda=k/h;
a=2;
u=zeros(nx,nt);
%icf
icf=@(x) 1 + sin(2*pi*x);
bc_xa=@(x)1;
bc_xb=@(x)1;
u(:,1)=icf(x);
u(1,1)=1;
u(nx,:)=bc_xb(t);

for j=2:nt
    for i=2:nx-1
        u(i,j)=u(i,j-1) - a*lambda*(u(i+1,j-1)-u(i-1,j-1))/2 + lambda*lambda*a*a*(u(i+1,j-1) - 2*u(i,j-1) + u(i-1,j-1))/2; ;
    end
    
    u(1,j)=(u(1,j-1) + u(2,j-1) + a*lambda*(u(2,j-1) - u(1,j-1)) - u(2,j)*(1+a*lambda))/(1-a*lambda);
end
[X,Y]=meshgrid(t,x);
figure
surf(X,Y,u);
xlabel('x');
ylabel('y');
zlabel('u');
title('(iv)Lax Wendroff Approx.');