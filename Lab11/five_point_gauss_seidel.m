function [U x y] = five_point_gauss_seidel(a_x,b_x,a_y,b_y,h,k,bc,F)
y = [a_y:k:b_y];
x = [a_x:h:b_x];
m = length(x);
n = length(y);
U = zeros(m,n);
lambda = h^2/k^2;
%Boundary Condition
j=1;
for i=1:m
    U(i,j) = bc(x(i),y(j));
end
j = n;
for i=1:m
    U(i,j) = bc(x(i),y(j));
end

i = 1;
for j=1:n
    U(i,j) = bc(x(i),y(j));
end
i = m;
for j=1:n
    U(i,j) = bc(x(i),y(j));
end
%%End of Boundary Condition
%Gauss - Seidel Iterations
tol = 1e-3;
curr_max_error = 100;
while curr_max_error > tol
    curr_max_error = 0;
    for i=2:m-1
        for j=2:n-1
            U_prev = U(i,j);
            U(i,j) = (U(i-1,j) + U(i+1,j) + lambda*U(i,j-1) + lambda*U(i,j+1))/(2+2*lambda) -(h^2)*F(x(i),y(j));
        curr_max_error = max(curr_max_error,abs(U_prev - U(i,j)));
        end
    end
end

end