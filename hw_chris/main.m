clear;
clc;

addpath('C:\Users\13695\Documents\Macro\compecon2011_64_20110718\compecon2011\CEtools')
%profile on
global beta r ro se gamma yPP ygrid
% note: using globals is bad practice! use structures instead


%declare parameters

beta=0.90;                         % Discounting
r=0.92/beta-1;                     % r(1+\beta)<1  
ro=0.70;                           % Coefficient of the AR(1) process
se=0.1;                            % SE of the AR(1) process
gamma=2;                           % Risk aversion

k=11;                              % nodes for the distribution of income shocks
[e,w]=rouwenhorst(ro,se,k);        % Rouwenhorst method
yPP=w;                             % Transition matrix we get 
ygrid=e';                          % Nodes of the income process

% Bounds for state space: we look for policy functions in the bound
ymin=min(e);                     % Upper bound of income process
ymax=max(e);                     % Lower bound of income process

amin = 0;                        % no borrowing
amax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it


% Declare function space to approximate a'(a,y)

n=[25,k];                        % Number of nodes in a space (25) and y space (k=11)

% Lower and higher bound for the state space (a,y)
smin=[amin,ymin];                % Lower bound of cartesian product of state space
smax=[amax,ymax];                % Upper bound of cartesian product of state space

scale=1/2;                       % Call the compecon, can change to 1/3 check for convergence 
                                 % simply to make the nodes denser close to
                                 % the kink
fspace=fundef({'spli',  nodeunif(n(1),(smin(1)-amin+.01).^scale,(smax(1)-amin+.01).^scale).^(1/scale)+amin-.01,0,3},...
              {'spli',ygrid,0,1});   % SPline - "nodeunif": Uniform to more dense grid        
% fspace is the guess of our policy function: a,y->x           

grid=funnode(fspace);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)
ns=length(s);

c=funfitxy(fspace,s,r/(1+r)*s(:,1)+exp(s(:,2)));                %guess that keep constant assets
% funfitxy:  Computes interpolation coefficients for d-dim function.

tic
for it=1:101                    % Run 100 iterations
cnew=c;                         % Temp value of coefficients 
%solve;      
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

%% Plot 
sfine=gridmake(nodeunif(n(1)*2,smin(1),smax(1)),ygrid);
xfine=funeval(c,fspace,sfine);
disp('mean and max errors')
disp(mean(abs(euler(xfine,c,fspace,sfine).*((1+r)*sfine(:,1)+exp(sfine(:,2))-xfine>amin+.01))))
disp(max(abs(euler(xfine,c,fspace,sfine).*((1+r)*sfine(:,1)+exp(sfine(:,2))-xfine>amin+.01))))

set(gca,'fontsize',8)
close all
figure(1)
subplot(2,1,1)
sfine=gridmake(nodeunif(n(1)*4,smin(1),smax(1)),0); %ygrid(floor(k/2)+2));
xfine=funeval(c,fspace,sfine);
plot(sfine(:,1),xfine)
xlabel({'$a$'},'Interpreter','latex')
ylabel({'$x(a,\bar{y})$'},'Interpreter','latex')
title({'Consumption policy function, $y=\bar{y}$'},'Interpreter','latex')
set(gca,'FontSize',8);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(n(2)*4,smin(2),smax(2)));
xfine=funeval(c,fspace,sfine);
plot(exp(sfine(:,2)),xfine)
xlabel('$e^{y}$','Interpreter','latex')
ylabel('$x(0,y)$','Interpreter','latex')
title({'Consumption policy function, $a=0$'},'Interpreter','latex')
set(gca,'FontSize',8);
print -djpeg -r600 hw_sample1


figure(2)
subplot(2,1,1)
sfine=gridmake(nodeunif(n(1)*4,smin(1),smax(1)),0);
xfine=funeval(c,fspace,sfine);
plot(sfine(:,1),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine)
xlabel('$a$','Interpreter','latex')
ylabel('$a^{\prime}(a,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $y=\bar{y}$'},'Interpreter','latex')
set(gca,'FontSize',8);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(n(2)*4,smin(2),smax(2)));
xfine=funeval(c,fspace,sfine);
plot(exp(sfine(:,2)),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine)
xlabel('$e^y$','Interpreter','latex')
ylabel('$a^{\prime}(0,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $a=0$'},'Interpreter','latex')
set(gca,'FontSize',8);
print -djpeg -r600 hw_sample2
