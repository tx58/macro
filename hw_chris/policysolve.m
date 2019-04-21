function c = policysolve(fspace, c)
% Input: number of grids
% Output: policy coefficient 
global beta r ro se gamma yPP ygrid wbar amin
% note: using globals is bad practice! use structures instead


%declare parameters     
grid=funnode(fspace);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)
ns=length(s);

tic
for it=1:301                    % Run 100 iterations
cnew=c;                         % Temp value of coefficients 
%solve;      s
    % Bisection method to solve the problem
    a=.01*ones(ns,1);
    b=(1+r)*s(:,1)+exp(s(:,2))-amin;  % Greatest level of consumption for each state
    tol=1e-8; %tolerance level

    fa=euler(a,c,fspace,s);       % Euler equation at the lower point
    fb=euler(b,c,fspace,s);       % Euler equation at the higher point

    x=zeros(ns,1);

    % Start bisection
    dx = 0.5*(b - a);

    x = a + dx;                       %  start consumption guess at midpoint
    sb=sign(fb);                      %  Sign at Higher point
    dx = sb.*dx;                      %  we increase or decrease x depending if f(b) is positive or negative                     

    i=0;
    
    % Inner loop: solve the nonlinear Euler equation for x
    while any(abs(dx)>tol)
       i=i+1;
        dx = 0.5*dx;
        x = x - sign(euler(x,c,fspace,s)).*dx;   % Here x is a vector, so we update the whole vector x each time 
    end

    x(fb>=0)=b(fb>=0);
    c=funfitxy(fspace,s,x);         % Find new coefficients given x
fprintf('%4i %6.2e\n',[it,norm(c-cnew)]);
if norm(c-cnew)<1e-7, break, end
end
toc  % Stop of time
%profile off
%profile viewer

end




