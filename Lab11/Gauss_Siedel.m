function y=Gauss_Siedel(A,b,n)
x=zeros(1,n);
tol=1e-2; itr=0;normVal=Inf;
while (normVal>tol)
    x_old=x; 
    for i=1:n
        sigma=0;
        for j=1:i-1
                sigma=sigma+A(i,j)*x(j);
        end
        for j=i+1:n
                sigma=sigma+A(i,j)*x_old(j);
        end   
        x(i)=(1/A(i,i))*(b(i)-sigma);
    end
    itr=itr+1;
    normVal=norm(x_old-x);
end
y=x;
y = y';
end