x = [0:6:84];
y = [124 134 148 156 147 133 121 109 99 85 78 89 104 116 123];
fprintf('Length of the track = %f feet\n',composite_trapezoid(x,y,length(x)));

function integral = composite_trapezoid(x,y,n)

integral = 0;
for i=1:n-1
    integral = integral + (x(i+1)-x(i))*(y(i)+y(i+1))/2;
end

end