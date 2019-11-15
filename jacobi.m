function output=jacobi(f,x,y,x0,x1,y0,y1)
    n=length(x)-1;
    h=x(2)-x(1);
    output=zeros(n+1,n+1);
    for i=1:n+1
        output(1,i)=y0(y(i));
        output(n+1,i)=y1(y(i));
        output(i,1)=x0(x(i));
        output(i,n+1)=x1(x(i));
    end
    temp=output;
    eps=10^-3;
    df=100;
    while df>eps
        df=0;
        for i=2:n
            for j=2:n
                temp(i,j)=0.25*(output(i-1,j)+output(i,j-1)+output(i,j+1)+output(i+1,j)-h*h*f(x(i),y(j)));
                df=max(df,abs(temp(i,j)-output(i,j)));
            end
        end
        output=temp;
    end
end