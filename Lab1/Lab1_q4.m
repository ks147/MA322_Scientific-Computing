a=1;
b=25;

%Q4(a) 
x=[];
i = 1;
x(1) = 1;
while(i < 5)
    x(i+1) = (1 + (7-x(i)^5)/(x(i)^2))^3;
    x(i+1) = x(i+1)*x(i);
    i = i+1;
end
i = 3;
rate = log(abs((x(i+1)-x(i))/(x(i)-x(i-1))));
rate = rate/(log(abs((x(i)-x(i-1))/(x(i-1)-x(i-2)))));
fprintf('Order of Convergence at 4th iteration for (a): %d \n',rate);

%Q4(b)
x=[];
i = 1;
x(1) = 1;
while(i < 5)
    x(i+1) = x(i) - (x(i)^5 - 7)/(x(i)^2);
    i = i+1;
end
i = 3;
rate = log(abs((x(i+1)-x(i))/(x(i)-x(i-1))));
rate = rate/(log(abs((x(i)-x(i-1))/(x(i-1)-x(i-2)))));
fprintf('Order of Convergence at 4th iteration for (b): %d \n',rate);

%Q4(c)
x=[];
i = 1;
x(1) = 1;
while(i < 5)
    x(i+1) = x(i) - (x(i)^5 - 7)/(5*(x(i)^4));
    i = i+1;
end
i = 3;
rate = log(abs((x(i+1)-x(i))/(x(i)-x(i-1))));
rate = rate/(log(abs((x(i)-x(i-1))/(x(i-1)-x(i-2)))));
fprintf('Order of Convergence at 4th iteration for (c): %d \n',rate);

%Q4(d)x=[];
i = 1;
x(1) = 1;
while(i < 5)
    x(i+1) = x(i) - (x(i)^5 - 7)/12;
    i = i+1;
end

i = 3;
rate = log(abs((x(i+1)-x(i))/(x(i)-x(i-1))));
rate = rate/(log(abs((x(i)-x(i-1))/(x(i-1)-x(i-2)))));
fprintf('Order of Convergence at 4th iteration for (c): %d \n',rate);

