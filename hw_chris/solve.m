% Bisection method to solve the problem
a=.01*ones(ns,1);
b=(1+r)*s(:,1)+exp(s(:,2))-amin;
tol=1e-8; %tolerance level

fa=euler(a,c,fspace,s,e,w);       % Euler equation at the lower point
fb=euler(b,c,fspace,s,e,w);       % Euler equation at the higher point

x=zeros(ns,1);

% Start bisection
dx = 0.5*(b - a);

x = a + dx;                       %  start at midpoint
sb=sign(fb);                      %  Sign at Higher point
dx = sb.*dx;                      %  we increase or decrease x depending if f(b) is positive or negative                     

i=0;
  while any(abs(dx)>tol)
   i=i+1;
    dx = 0.5*dx;
    x = x - sign(euler(x,c,fspace,s,e,w)).*dx;   
  end


x(fb>=0)=b(fb>=0);
