x = [1 5/4 3/2 7/4 2];
y = [10 8 7 6 5];
fprintf('Area under f(x) using composite Trapezoid = %f\n',composite_trapezoid(x,y,length(x)));

function integral = composite_trapezoid(x,y,n)

integral = 0;
for i=1:n-1
    integral = integral + (x(i+1)-x(i))*(y(i)+y(i+1))/2;
end

end