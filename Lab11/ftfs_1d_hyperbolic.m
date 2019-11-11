function [U x t] = ftfs_1d_hyperbolic(a,a_x,b_x,a_t,b_t,h,k,ic,bc_2)
R = a*k/h;

x = [a_x:h:b_x];    %Discretize space
t = [a_t:k:b_t];    %Discretize time range
n = length(t);
m = length(x);
U = zeros(m,n);
%Initial Condition
U(:,1) = ic(x);

%Boundary Condition
U(m-1,:) = bc_2(t);

%Only for 3 time levels
for j=2:3    
    for i=2:m-1
        U(i,j) = (1+R)*U(i,j-1) + R*U(i+1,j);
    end
end

end