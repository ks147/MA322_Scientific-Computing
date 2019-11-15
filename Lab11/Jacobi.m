function y=Jacobi(A,b,n)
x=zeros(1,n);
tol=1e-3; itr=0;normVal=Inf;
while normVal>tol
    xold=x;
    for i=1:n
        sigma=0;  
        for j=1:n          
            if j~=i
                sigma=sigma+A(i,j)*x(j);
            end          
        end
        x(i)=(1/A(i,i))*(b(i)-sigma);
    end
    
    itr=itr+1;
    normVal=abs(xold-x);
end
y=x; 
end