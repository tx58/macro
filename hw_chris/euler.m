function fval = euler(x,c,fspace,s)

global beta r ro gamma yPP ygrid

ns=size(s,1); % Length of the state space
a=s(:,1); % Current period asset: state
y=exp(s(:,2)); % Current y: state
aprime=(1+r)*a+y-x; % Find next period asset from Budget constraint
fval=x.^(-gamma); % u' function
as=ns/length(ygrid); 

for i=1:length(ygrid)
  xprime=funeval(c,fspace,[aprime,ygrid(i)*ones(ns,1)]); % Get next period consumption, Evaluate functions at linear bases
  for j=1:length(ygrid)
    ypp(1+(j-1)*as:j*as,1)=yPP(j,i); % Transition matrix
  end 
    fval=fval-beta*(1+r)*ypp.*xprime.^(-gamma); % Euler equation
end
