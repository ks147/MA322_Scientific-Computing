function [U x y] = five_point_scheme(a_x,b_x,a_y,b_y,h,k,bc,F,part,method)

x = [a_x:h:b_x];
y = [a_y:k:b_y];
m = length(x);
n = length(y);
lambda = h^2/k^2;
A = zeros(m*n,m*n);
%lth row in A is the equation for 
%U(i,j) s.t. l = i + (j-1)m
%i--columns,x
%j--rows,y
%U(x,y)
%Boundary Condition
f = zeros(m*n,1);
j = 1;
for i=1:m
    A(i+(j-1)*m,i+(j-1)*m) = 1;
    f(i+(j-1)*m) = bc(x(i),y(j));
end
j = n;
for i=1:m
    A(i+(j-1)*m,i+(j-1)*m) = 1;
    f(i+(j-1)*m) = bc(x(i),y(j));
end
i = 1;
for j=1:n
    if part=='c'
        A(i+(j-1)*m,i+1+(j-1)*m) = 1/h;
        A(i+(j-1)*m,i+(j-1)*m) = 1 - 1/h;
        f(i+(j-1)*m) = 2 - y(j);
    else
        A(i+(j-1)*m,i+(j-1)*m) = 1;
        f(i+(j-1)*m) = bc(x(i),y(j));
    end
end
i = m;
for j=1:n
    A(i+(j-1)*m,i+(j-1)*m) = 1;
    f(i+(j-1)*m) = bc(x(i),y(j));
end
%end of BC
%five point scheme
for i=2:m-1
    for j=2:n-1
       %-4u(i,j) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) = h^2*F(x(i),y(j))
       if part == 'd'
       A(i+(j-1)*m,i+(j-1)*m) = h^2 - 2*h -4;
       A(i+(j-1)*m,i+(j-2)*m) = 1;
       A(i+(j-1)*m,i+(j)*m) = 1+h;
       A(i+(j-1)*m,i-1+(j-1)*m) = 1;
       A(i+(j-1)*m,i+1+(j-1)*m) = 1+h;
       else
       A(i+(j-1)*m,i+(j-1)*m) = -2*(lambda + 1);
       A(i+(j-1)*m,i+(j-2)*m) = lambda;
       A(i+(j-1)*m,i+(j)*m) = lambda;
       A(i+(j-1)*m,i-1+(j-1)*m) = 1;
       A(i+(j-1)*m,i+1+(j-1)*m) = 1;
       end
       f(i+(j-1)*m) = (h.^2)*F(x(i),y(j));
    end
end
if method=='Gauss Elimination'
    W = A\f;
elseif method=='Gauss Seidel'
    W = Gauss_Siedel(A,f,m*n);
else 
    W = Jacobi(A,f,m*n);
end

U = reshape(W,n,m);

% for i = 1:m
%     for j = 1:n
%         U(i,j) = W(i+(j-1)*m);
%     end
% end

end
