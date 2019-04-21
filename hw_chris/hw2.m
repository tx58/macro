% Tianli Xia
% Macro Chris Homework 2
% Dynamic programming with kink in Value function
clear 
addpath('C:\Users\13695\Documents\Macro\compecon2011_64_20110718\compecon2011\CEtools')
global beta r ro se gamma yPP ygrid wbar amin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Q1: Tacuchen and Rouwnhorst Method
% w_t= w_bar + \rho*w_t-1 + \eps_t
% w_t will follow a N(w_bar/(1-rho), \sigma2_eps/(1-rho^2))
% y_t follows Lognormal distribution with mean 1: exp(\mu + \sigma^{2}/2) =1
se2= 0.06; 
se= sqrt(se2);             % SE of the AR(1) process
ro= 0.9;                   % Coefficient of the AR(1) process
wbar= -se2/2/(1+ ro);

% 1. Tauchen Method
m= 11; % number of discrete points we approximate
sigmaV= sqrt( se2/(1-ro^2) );
fprintf("True se: \n %f \n", sigmaV);

P= ones(m);
d= 6*sigmaV/(m-1);
v= zeros(m,1);
v= -3*sigmaV : d : 3*sigmaV;
w= v+ wbar/(1-ro);

for j=1:m % Note that in the loop i=1-> Y(i)=-inf, so we need P(1,.)=Y(2)
    for k=1:m
        if k==1
            P(j,k)=normcdf( (v(k)+d/2-ro*v(j)) /se );
        elseif k==m
            P(j,k)=1-normcdf( (v(k)-d/2-ro*v(j)) /se );
        else
            P(j,k)=normcdf( (v(k)+d/2-ro*v(j)) /se )-...
               normcdf( (v(k)-d/2-ro*v(j)) /se );
        end
    end
end

K=10000;
what(1)=3;

[V lamda]=eigs(P');  % Find eigenvectors, solve for PS=S.
S=V(:,1)';
S=S/sum(S); % Solve for the steady state of the transition matrix
fprintf("The mean of Tauchen method is:\n %f \n", sum(exp(w).*S));
% var(exp(w), (S))
fprintf("The se of Tauchen method is:\n %f \n", sqrt(var(w, S)));

what(1)=(m+1)/2;
for i=1:K-1 
     what(i+1)= 1+sum(rand(1)>cumsum(P(what(i),:)));
end
fprintf("The ro of Tauchen method is:\n %f \n", corr(exp(w(what(2:end)))',exp(w(what(1:end-1)))'));


% for i=1:K-1 
%     what(i+1)= 1+sum(rand(1)>cumsum(P(what(i),:)));
% end
% mdl = arima(1,0,0); % 2 the lag order
% EstMdl = estimate(mdl,w(what(1001:K))'); % y is your data

%% 2. rouwenhorst
m=11;
[zgrid, P2] = rouwenhorst(ro, se, m);
wgrid= zgrid+ wbar/(1-ro);
K=10000;
[V2 ~]=eig(P2');  % Find eigenvectors, solve for PS=S.
S2=V2(:,1)';
S2=S2/sum(S2); % Solve for the steady state of the transition matrix
fprintf("The mean of rouwenhorst method is:\n %f \n", sum(exp(wgrid).*S2));
% var(exp(w), (S))
fprintf("The se of rouwenhorst method is:\n %f \n", sqrt(var(wgrid, S2)));
what(1)= (m+1)/2;
for i=1:K-1 
     what(i+1)= 1+sum(rand(1)>cumsum(P(what(i),:)));
end
fprintf("The ro of Tauchen method is:\n %f \n", corr(exp(wgrid(what(2:end)))',exp(wgrid(what(1:end-1)))'));

% mdl = arima(1,0,0); % 2 the lag order
% EstMdl = estimate(mdl,wgrid(what(1001:K))'); % y is your data

%% 3. Policy Function
%declare parameters
l=25;
k=9;
beta=0.95;                         % Discounting
r=0.02;                            % r(1+\beta)<1   
gamma=1;                           % Risk aversion
se2= 0.06; 
se= sqrt(se2);             % SE of the AR(1) process
ro= 0.9;                   % Coefficient of the AR(1) process
wbar= -se2/2/(1+ ro);
[e,w]=rouwenhorst(ro, se, k);  %  Rouwenhorst method
yPP=w;                         % Transition matrix we get 
ygrid=(e'+wbar/(1-ro));        % Nodes of the income process

% Bounds for state space: we look for policy functions in the bound
ymin=min(ygrid);                     % Upper bound of income process
ymax=max(ygrid);                     % Lower bound of income process

amin = 0;                        % no borrowing
amax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it

% Declare function space to approximate a'(a,y)
n=[l,k];                        % Number of nodes in a space (25) and y space (k=11)

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
c=funfitxy(fspace,s,r/(1+r)*s(:,1)+exp(s(:,2)));                %guess that keep constant assets
% funfitxy:  Computes interpolation coefficients for d-dim function.

c1  = policysolve(fspace, c);
gamma=2; 
c2  = policysolve(fspace, c);                       
gamma=3;                           
c3  = policysolve(fspace, c);

%% Compare in the graph
figure(1)
subplot(2,1,1)
sfine=gridmake(nodeunif(n(1)*4,smin(1),smax(1)),0); %ygrid(floor(k/2)+2));
xfine1=funeval(c1,fspace,sfine);
xfine2=funeval(c2,fspace,sfine);
xfine3=funeval(c3,fspace,sfine);
plot((sfine(:,1)),xfine1)
hold on
plot((sfine(:,1)),xfine2)
hold on
plot((sfine(:,1)),xfine3)
xlabel({'$a$'},'Interpreter','latex')
ylabel({'$x(a,\bar{y})$'},'Interpreter','latex')
title({'Consumption policy function, $y=\bar{y}$'},'Interpreter','latex')
legend("\gamma=1","\gamma=2","\gamma=3")
set(gca,'FontSize',12);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(n(2)*4,smin(2),smax(2)));
xfine1=funeval(c1,fspace,sfine);
xfine2=funeval(c2,fspace,sfine);
xfine3=funeval(c3,fspace,sfine);
plot(exp(sfine(:,2)),xfine1)
hold on
plot(exp(sfine(:,2)),xfine2)
hold on
plot(exp(sfine(:,2)),xfine3)
xlabel('$e^{y}$','Interpreter','latex')
ylabel('$x(0,y)$','Interpreter','latex')
title({'Consumption policy function, $a=0$'},'Interpreter','latex')
legend("\gamma=1","\gamma=2","\gamma=3")
set(gca,'FontSize',12);
print -djpeg -r600 hw_gamma_consumption

figure(2)
subplot(2,1,1)
sfine=gridmake(nodeunif(n(1)*4,smin(1),smax(1)),0);
xfine1=funeval(c1,fspace,sfine);
xfine2=funeval(c2,fspace,sfine);
xfine3=funeval(c3,fspace,sfine);
plot(sfine(:,1),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine1)
hold on
plot(sfine(:,1),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine2)
hold on
plot(sfine(:,1),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine3)
legend("\gamma=1","\gamma=2","\gamma=3")
xlabel('$a$','Interpreter','latex')
ylabel('$a^{\prime}(a,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $y=\bar{y}$'},'Interpreter','latex')
set(gca,'FontSize',12);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(n(2)*4,smin(2),smax(2)));
xfine1=funeval(c1,fspace,sfine);
xfine2=funeval(c2,fspace,sfine);
xfine3=funeval(c3,fspace,sfine);
plot(exp(sfine(:,2)),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine1)
hold on
plot(exp(sfine(:,2)),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine2)
hold on
plot(exp(sfine(:,2)),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine3)
legend("\gamma=1","\gamma=2","\gamma=3")
xlabel('$e^y$','Interpreter','latex')
ylabel('$a^{\prime}(0,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $a=0$'},'Interpreter','latex')
set(gca,'FontSize',12);
print -djpeg -r600 hw_gamma_saving
%% Simulate a path of income shocks using the transition matrix
% Generate several shocks
w= zeros(K,1);
w(1)= (k+1)/2;
K=1000;
for t=1:K-1
    w(t+1)= 1+sum(rand(1)>cumsum(yPP(w(t),:)));
end
[w1, a1, x1]= simulate(fspace, c1, K, w);
[w2, a2, x2]= simulate(fspace, c2, K, w);
[w3, a3, x3]= simulate(fspace, c3, K, w);
period=1:1:1000;
plot(period, x1)
hold on
plot(period, x2)
hold on
plot(period, x3)
legend("\gamma=1","\gamma=2","\gamma=3")
xlabel('Periods','Interpreter','latex')
ylabel('Consumption','Interpreter','latex')
title({'Consumption with simulated path'},'Interpreter','latex')
print -djpeg -r600 hw_gamma_simulation

str="Gamma Stats";
sheet=1;
A=[1                1                 1
  var(x1)           var(x2)           var(x3)
  corr(ygrid(w),x1) corr(ygrid(w),x2) corr(ygrid(w),x3)];
xlswrite(str, A, sheet, 'B2' );

%% Vary across se
% se2=0.01
gamma=2;
se2= 0.01; 
se= sqrt(se2);             % SE of the AR(1) process
ro= 0.9;                   % Coefficient of the AR(1) process
wbar= -se2/2/(1+ ro);
[e,w]=rouwenhorst(ro, se, k);  %  Rouwenhorst method
yPP=w;                         % Transition matrix we get 
ygrid=(e'+wbar/(1-ro));        % Nodes of the income process

% Bounds for state space: we look for policy functions in the bound
ymin=min(ygrid);                     % Upper bound of income process
ymax=max(ygrid);                     % Lower bound of income process

amin = 0;                        % no borrowing
amax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it

% Declare function space to approximate a'(a,y)
n=[l,k];                        % Number of nodes in a space (25) and y space (k=11)

% Lower and higher bound for the state space (a,y)
smin=[amin,ymin];                % Lower bound of cartesian product of state space
smax=[amax,ymax];                % Upper bound of cartesian product of state space

scale=1/2;                       % Call the compecon, can change to 1/3 check for convergence 
                                 % simply to make the nodes denser close to
                                 % the kink
fspace1=fundef({'spli',  nodeunif(n(1),(smin(1)-amin+.01).^scale,(smax(1)-amin+.01).^scale).^(1/scale)+amin-.01,0,3},...
              {'spli',ygrid,0,1});   % SPline - "nodeunif": Uniform to more dense grid        
% fspace is the guess of our policy function: a,y->x           

grid=funnode(fspace1);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)
c=funfitxy(fspace1,s,r/(1+r)*s(:,1)+exp(s(:,2)));                %guess that keep constant assets
% funfitxy:  Computes interpolation coefficients for d-dim function.
c4  = policysolve(fspace1, c);


% se2=0.12
se2= 0.12; 
se= sqrt(se2);             % SE of the AR(1) process
ro= 0.9;                   % Coefficient of the AR(1) process
wbar= -se2/2/(1+ ro);
[e,w]=rouwenhorst(ro, se, k);  %  Rouwenhorst method
yPP=w;                         % Transition matrix we get 
ygrid=(e'+wbar/(1-ro));        % Nodes of the income process

% Bounds for state space: we look for policy functions in the bound
ymin=min(ygrid);                     % Upper bound of income process
ymax=max(ygrid);                     % Lower bound of income process

amin = 0;                        % no borrowing
amax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it

% Declare function space to approximate a'(a,y)
n=[l,k];                        % Number of nodes in a space (25) and y space (k=11)

% Lower and higher bound for the state space (a,y)
smin=[amin,ymin];                % Lower bound of cartesian product of state space
smax=[amax,ymax];                % Upper bound of cartesian product of state space

scale=1/2;                       % Call the compecon, can change to 1/3 check for convergence 
                                 % simply to make the nodes denser close to
                                 % the kink
fspace2=fundef({'spli',  nodeunif(n(1),(smin(1)-amin+.01).^scale,(smax(1)-amin+.01).^scale).^(1/scale)+amin-.01,0,3},...
              {'spli',ygrid,0,1});   % SPline - "nodeunif": Uniform to more dense grid        
% fspace is the guess of our policy function: a,y->x           

grid=funnode(fspace2);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)
c=funfitxy(fspace2,s,r/(1+r)*s(:,1)+exp(s(:,2)));                %guess that keep constant assets
% funfitxy:  Computes interpolation coefficients for d-dim function.

c5  = policysolve(fspace2, c);

%% Plot
figure(2)
subplot(2,1,1)
sfine=gridmake(nodeunif(n(1)*3,smin(1),smax(1)),0);
xfine1=funeval(c2,fspace,sfine);
xfine4=funeval(c4,fspace1,sfine);
xfine5=funeval(c5,fspace2,sfine);
plot(sfine(:,1),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine4)
hold on
plot(sfine(:,1),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine1)
hold on
plot(sfine(:,1),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine5)
legend("\sigma^2_\epsilon = 0.01", "\sigma^2_\epsilon = 0.06","\sigma^2_\epsilon = 0.12","location", "b")
xlabel('$a$','Interpreter','latex')
ylabel('$a^{\prime}(a,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $y=\bar{y}$'},'Interpreter','latex')
set(gca,'FontSize',12);
subplot(2,1,2)
sfine=gridmake(0,nodeunif(n(2)*3,-1,1));
xfine1=funeval(c2,fspace,sfine);
xfine4=funeval(c4,fspace1,sfine);
xfine5=funeval(c5,fspace2,sfine);
plot(exp(sfine(:,2)),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine4)
hold on
plot(exp(sfine(:,2)),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine1)
hold on
plot(exp(sfine(:,2)),(1+r)*sfine(:,1)+exp(sfine(:,2))-xfine5)
legend("\sigma^2_\epsilon = 0.01", "\sigma^2_\epsilon = 0.06","\sigma^2_\epsilon = 0.12","location", "b")
xlabel('$e^y$','Interpreter','latex')
ylabel('$a^{\prime}(0,\bar{y})$','Interpreter','latex')
title({'Savings policy function, $a=0$'},'Interpreter','latex')
set(gca,'FontSize',12);
print -djpeg -r600 hw_se_saving


%% Natural debt limit: no borrowing constraint

l=25;
k=9;
beta=0.95;                         % Discounting
r=0.02;                            % r(1+\beta)<1   
gamma=2;                           % Risk aversion
se2= 0.06; 
se= sqrt(se2);             % SE of the AR(1) process
ro= 0.9;                   % Coefficient of the AR(1) process
wbar= -se2/2/(1+ ro);
[e,w]=rouwenhorst(ro, se, k);  %  Rouwenhorst method
yPP=w;                         % Transition matrix we get 
ygrid=(e'+wbar/(1-ro));        % Nodes of the income process

% Bounds for state space: we look for policy functions in the bound
ymin=min(ygrid);                     % Upper bound of income process
ymax=max(ygrid);                     % Lower bound of income process

amin = -exp(ygrid(2))/r;                        % no borrowing
amax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it

% Declare function space to approximate a'(a,y)
n=[l,k];                        % Number of nodes in a space (25) and y space (k=11)

% Lower and higher bound for the state space (a,y)
smin=[amin,ymin];                % Lower bound of cartesian product of state space
smax=[amax,ymax];                % Upper bound of cartesian product of state space

scale=1/2;                       % Call the compecon, can change to 1/3 check for convergence 
                                 % simply to make the nodes denser close to
                                 % the kink
fspace3=fundef({'spli',  nodeunif(n(1),(smin(1)-amin+.01).^scale,(smax(1)-amin+.01).^scale).^(1/scale)+amin-.01,0,3},...
              {'spli',ygrid,0,1});   % SPline - "nodeunif": Uniform to more dense grid        
% fspace is the guess of our policy function: a,y->x           

grid=funnode(fspace3);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)
c=funfitxy(fspace3,s,r/(1+r)*s(:,1)+exp(s(:,2)));                %guess that keep constant assets
% funfitxy:  Computes interpolation coefficients for d-dim function.

c6  = policysolve(fspace3, c);
%% Plot
figure(1)
subplot(2,1,1)
sfine=gridmake(nodeunif(n(1)*4,0,smax(1)),0); %ygrid(floor(k/2)+2));
xfine2=funeval(c2,fspace,sfine);
xfine6=funeval(c6,fspace3,sfine);
plot((sfine(:,1)),xfine2)
hold on
plot((sfine(:,1)),xfine6)
xlabel({'$a$'},'Interpreter','latex')
ylabel({'$x(a,\bar{y})$'},'Interpreter','latex')
title({'Consumption policy function, $y=\bar{y}$'},'Interpreter','latex')
legend("With zero constraint", "With nature constraint", "location", "b")
set(gca,'FontSize',12);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(n(2)*4,smin(2),smax(2)));
xfine2=funeval(c2,fspace,sfine);
xfine6=funeval(c6,fspace3,sfine);
plot(exp(sfine(:,2)),xfine2)
hold on
plot(exp(sfine(:,2)),xfine6)
xlabel('$e^{y}$','Interpreter','latex')
ylabel('$x(0,y)$','Interpreter','latex')
title({'Consumption policy function, $a=0$'},'Interpreter','latex')
legend("With zero constraint", "With nature constraint", "location", "b")
set(gca,'FontSize',12);
print -djpeg -r600 hw_borrowing_constraint_consumption


%% Compare with no borrowing constraints

l=25;
k=9;
beta=0.95;                         % Discounting
r=0.02;                            % r(1+\beta)<1   
gamma=2;                           % Risk aversion
se2= 0.06; 
se= sqrt(se2);             % SE of the AR(1) process
ro= 0.9;                   % Coefficient of the AR(1) process
wbar= -se2/2/(1+ ro);
[e,w]=rouwenhorst(ro, se, k);  %  Rouwenhorst method
yPP=w;                         % Transition matrix we get 
ygrid=(e'+wbar/(1-ro));        % Nodes of the income process

% Bounds for state space: we look for policy functions in the bound
ymin=min(ygrid);                     % Upper bound of income process
ymax=max(ygrid);                     % Lower bound of income process

amin = -10*exp(ymin);                        % no borrowing
amax = 10*exp(ymax);             % guess an upper bound on a, check later that do not exceed it

% Declare function space to approximate a'(a,y)
n=[l,k];                        % Number of nodes in a space (25) and y space (k=11)

% Lower and higher bound for the state space (a,y)
smin=[amin,ymin];                % Lower bound of cartesian product of state space
smax=[amax,ymax];                % Upper bound of cartesian product of state space

scale=1/3;                       % Call the compecon, can change to 1/3 check for convergence 
                                 % simply to make the nodes denser close to
                                 % the kink
fspace4=fundef({'spli',  nodeunif(n(1),(smin(1)-amin+.01).^scale,(smax(1)-amin+.01).^scale).^(1/scale)+amin-.01,0,3},...
              {'spli',ygrid,0,1});   % SPline - "nodeunif": Uniform to more dense grid        
% fspace is the guess of our policy function: a,y->x           

grid=funnode(fspace4);
s=gridmake(grid); %collection of  states (all a with y1... all a with y2... and so on)
c=funfitxy(fspace4,s,r/(1+r)*s(:,1)+exp(s(:,2)));                %guess that keep constant assets
% funfitxy:  Computes interpolation coefficients for d-dim function.

c7  = policysolve(fspace4, c);


%% Plot and find insurance level 
figure(1)
subplot(2,1,1)
sfine=gridmake(nodeunif(n(1)*4,0,smax(1)),0); %ygrid(floor(k/2)+2));
xfine7=funeval(c7,fspace4,sfine);
xfine6=funeval(c6,fspace3,sfine);
plot((sfine(:,1)),xfine7)
hold on
plot((sfine(:,1)),xfine6)
xlabel({'$a$'},'Interpreter','latex')
ylabel({'$x(a,\bar{y})$'},'Interpreter','latex')
title({'Consumption policy function, $y=\bar{y}$'},'Interpreter','latex')
legend("With nature constraint", "No constraint", "location", "b")
set(gca,'FontSize',12);

subplot(2,1,2)
sfine=gridmake(0,nodeunif(n(2)*4,smin(2),smax(2)));
xfine7=funeval(c7,fspace4,sfine);
xfine6=funeval(c6,fspace3,sfine);
plot(exp(sfine(:,2)),xfine7)
hold on
plot(exp(sfine(:,2)),xfine6)
xlabel('$e^{y}$','Interpreter','latex')
ylabel('$x(0,y)$','Interpreter','latex')
title({'Consumption policy function, $a=0$'},'Interpreter','latex')
legend("With nature constraint", "No constraint", "location", "b")
set(gca,'FontSize',12);
print -djpeg -r600 hw_no_borrowing_constraint_consumption

w= zeros(K,1);
w(1)= (k+1)/2;
K=1000;
eps(1)=0;
for t=1:K-1
    w(t+1)= 1+sum(rand(1)>cumsum(yPP(w(t),:)));
    eps(t+1)= ygrid(w(t+1))-0.9* ygrid(w(t));
end
[w6, a6, x6]= simulate(fspace3, c6, K, w);
[w7, a7, x7]= simulate(fspace4, c6, K, w);

insurance6= 1- cov(eps, log(x6))/ var(eps)
insurance7= 1- cov(eps, log(x7))/ var(eps)

xlswrite(str, [insurance6(1,2) insurance7(1,2)], 3, 'B2' );

%% Cash-in-hand
% Simulate a normal distribution
[x,w] = qnwnorm(7, -(0.06)/2, (0.06));
sum(exp(x).*w)

