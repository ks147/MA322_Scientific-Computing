clear all
clc

a=6;
b=0;
m=11;
for j=1:11
  subplot(4,3,j)  
  x(1,j)=j-1;
  for i=1:10
    x(i+1,j)= mod(a*x(i,j)+b,m);
    u(i+1,j)= (x(i+1,j)/m);
  end;
  hist(u,11)
  disp(u)
end;
a=3;
for j=1:11
  subplot(4,3,j)  
  x(1,j)=j-1;
  for i=1:10
    x(i+1,j)= mod(a*x(i,j)+b,m);
    u(i+1,j)= (x(i+1,j)/m);
  end;
  hist(u,11)
  disp(u)
end;