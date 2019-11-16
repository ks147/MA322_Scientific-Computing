f = @(t,y) -2.*y;
exact_sol = @(t) exp(2.*t);
a = 0;
b = 1;
h = 1e-3;
alpha = 1;
milne(a,b,h,alpha,f,exact_sol)