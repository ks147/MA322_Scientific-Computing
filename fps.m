function output=fps(f,x,y,x0,x1,y0,y1)
    n=length(x)-1;
    h=x(2)-x(1);
    m=(n+1).^2;
    S=zeros(m,m);
    b=zeros(m,1);
    output=zeros(n+1,n+1);
    for i=1:n+1
        S(i,i)=1;
        S(m-i+1,m-i+1)=1;
    end
    v=zeros(m,1);
    for i=n+2:m-(n+1)
        if (mod(i,n+1)==1)
            S(i,i)=1;
            j=(i-1)/(n+1);
            b(i)=y0(y(n+1-j));
            continue;
        end
        if(mod(i,n+1)==0)
            S(i,i)=1;
            j=(i/(n+1))-1;
            b(i)=y1(y(n+1-j));
            continue;
        end
        S(i,i)=4;
        S(i,i+1)=-1;
        S(i,i-1)=-1;
        S(i,i+n+1)=-1;
        S(i,i-n-1)=-1;
        z=mod(i,n+1);
        xx=z;
        yy=(n+1)-(i-z)/(n+1);
        b(i)=h*h*f(x(xx),y(yy));
    end

    for i=1:(n+1)
        b(i)=x1(x(i));
        b(m-i+1)=x0(x(n+2-i));
    end
    v=S\b;
    for i=1:m
        z=mod(i,n+1);
        if z==0
            output((n+2)-((i-z)/(n+1)),n+1)=v(i);
            continue;
        end
        output((n+1)-((i-z)/(n+1)),z)=v(i);
    end
end