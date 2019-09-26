clear all
clc

m=244944;
a=1597;
b=0;
for j=1:5
    subplot(2,3,j);
    x(1,j)=10*j;
    for i=1:m-1
        x(i+1,j)=mod(a*x(i,j) +b, m);
        u(i+1,j)=x(i+1,j)/m;
    end;
    hist(j)=histogram(u(:,j),[0:0.05:1]);
end;

for i=1:5
    freq(:,i) = hist(i).Values;
end;
 disp(freq)
        