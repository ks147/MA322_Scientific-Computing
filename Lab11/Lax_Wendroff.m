function [U x t] = Lax_Wendroff(a_x,b_x,a_t,b_t,h,k,a,ic,bc)
%Ut = -aUx
x = [a_x:h:b_x];
t = [a_t:k:b_t];
m = length(x);
n = length(t);
U = zeros(m,n);  %U(x,t);
U(:,1) = ic(x);
U(m,:) = bc(t);
U(1,:) = bc(t);
R = a*k/h;

for j = 2:n
    for i=2:m-1
        U(i,j) = U(i,j-1) - 0.5*R*(U(i+1,j-1) - U(i-1,j-1)) + (R*R*0.5)*(U(i+1,j-1) - 2*U(i,j-1) + U(i-1,j-1));
    end
end


end
