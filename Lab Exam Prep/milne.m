function milne(a,b,h,alpha,f,exact_sol)
t = [a:h:b];
n = length(t);
y = zeros(n,1);
y(1) = alpha;
y(2) = y(1) + h*f(t(1),y(1));
w_bash = [3/2 -1/2];
w_milne = [1/3 4/3 1/3];
for i=3:n
    F = [];
    for j=1:2
        F(j) = f(t(i-j),y(i-j));
    end
    y(i) = y(i-1) + h*dot(w_bash,F); %Bashforth prediction
    for j=0:2
        F(j+1) = f(t(i-j),y(i-j));
    end
    y(i) = y(i-1) + h*dot(w_milne,F);    
end
error = abs(exact_sol(t) - y');
figure;
plot(t,error);

end