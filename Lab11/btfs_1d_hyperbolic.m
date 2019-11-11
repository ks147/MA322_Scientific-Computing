function [U x t] = btfs_1d_hyperbolic(a,a_x,b_x,a_t,b_t,h,k,ic,bc_2)
R = a*k/h;

x = [a_x:h:b_x];    %Discretize space
t = [a_t:k:b_t];    %Discretize time range
n = length(t);
m = length(x);
U = zeros(m,n);
U(:,1) = ic(x);
%IMPLICIT SCHEME
A = zeros(m,m);
%A(1,1) = 1;
A(m,m) = 1;

%Boundary Condition
U(m-1,:) = bc_2(t);

for i = 1:m-1
    A(i,i) = (1-R);
    A(i,i+1) = R;
end

%Solving only for 3 time levels 
%for solution at all time levels replace 3 by n
for j=2:3
    U(:,j) = A\U(:,j-1);
end

end