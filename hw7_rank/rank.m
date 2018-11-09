%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HW7 New Keynsian Model
% Tianli Xia
% November 5th, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializtion
clear
sigma= 2; % consumption elasiticity 
phi= 3; % labor elasticity
epsilon= 5; % different goods bundle elasticity
alpha= 0.3; % production coefficient
rhoa= 0.8; % technology AR(1) coefficient
sigmaa2= 1; % technology AR(1) residuals
beta= 0.99; % intertemporal discount factor
rho= -log(beta); % linearied beta
theta= 0.75; % ratio of firms which cannot change prices
psipi= 1.5; % taylor rule of inflation
psiy= 0.5; % taylor rule of output
T=100; % total period

% Some definitions:
psiya= (1+phi)/(sigma*(1-alpha)+phi+alpha); % tech effect on output
lamda= (1-theta)*(1-beta*theta)/theta*(1-alpha)/(1-alpha+alpha*epsilon);
kappa= lamda*(sigma+(phi+alpha)/(1-alpha));
gamma= 1/((1-beta*rhoa)*(sigma*(1-rhoa)+psiy)+ kappa*(psipi-rhoa));




%% A one-time Shock
% 1. Specify the shock
v= zeros(T);
v(2)= 1;
a(1)= 0; % initial technology level

% 2. Calculate the deviation of output and inflation
output(1)= 0; % start from steady state
inflation(1)=0;
ir(1)=rho;
r(1)= ir(1);

% 3. price-related variables:
price(1)= theta*inflation(1)/(1-theta); % price deviation from optimal p1*-pt
product(1)= -epsilon*(price(1))+ output(1); % output for product j: y_{t}(j)=-?(p_{t}^{j}-p_{t+k})+y_{t+k}
mc(1)= (sigma+(phi+alpha)/(1-alpha))*product(1); % mc for product j: mc_{t}(j)-mc_{t}= (?+((?+?)/(1-?)))(y_{t}(j)-y_{t})

% 4. iterate to get future response
for i=2:T
    a(i)= rhoa*a(i-1)+ v(i);
    output(i)= rhoa*output(i-1) - sigma* psiya* (1-rhoa)* (1-beta*rhoa)*gamma*v(i);
    inflation(i)= rhoa*inflation(i-1) - sigma* psiya* (1-rhoa)* kappa *gamma*v(i);
    ir(i)=rho+ psipi*inflation(i) + psiy* output(i);
    y(i)=psiya*a(i)+ output(i);
    r(i)=ir(i)-inflation(i);
    price(i)= price(i-1)- output(i);
    product(i)= -epsilon*(price(i))+ output(i);
    mc(i)=  (sigma+(phi+alpha)/(1-alpha))*product(i);
end

% 3. Consider the level too
outputn= psiya.*a+ (1-alpha)*(0-log(1-alpha))/((1-alpha)*sigma+ phi+ alpha);
inflationn= zeros(T);
irn= rho-sigma*(1-rhoa)*psiya.*a;

% 4. Plot response to a one-time shock
t=20;
x =2:t;
subplot(3,3,1);
plot(x,a(2:t))
title('technology')
subplot(3,3,2);
plot(x,output(2:t))
title('output gap')
subplot(3,3,3);
plot(x,inflation(2:t))
title('inflation')
subplot(3,3,4);
plot(x,ir(2:t))
title('nominal rate')
subplot(3,3,5);
plot(x,r(2:t))
title('real rate')
subplot(3,3,6);
plot(x,y(2:t))
title('output')
subplot(3,3,7);
plot(x,price(2:t))
title('price deviation')
subplot(3,3,8);
plot(x,product(2:t))
title('product(j)')
subplot(3,3,9);
plot(x,mc(2:t))
title('marginal cost(j)')
suptitle('Impulse Response for a positive technology shock')
print -djpeg -r600 q4_1.jpg

%% Plot 100 period simulation of the economy
% 1. Generate randomness technology shock:
v= sqrt(sigmaa2)*randn(100); % tech shock
a(1)= 0; % initial technology level

% 2. Calculate the deviation of output and inflation
output(1)= 0; % start from steady state
inflation(1)=0;
ir(1)=rho;
r(1)= ir(1);

% 3. price-related variables:
price(1)= theta*inflation(1)/(1-theta); % price deviation from optimal p1*-pt
product(1)= -epsilon*(price(1))+ output(1); % output for product j: y_{t}(j)=-?(p_{t}^{j}-p_{t+k})+y_{t+k}
mc(1)= (sigma+(phi+alpha)/(1-alpha))*product(1); % mc for product j: mc_{t}(j)-mc_{t}= (?+((?+?)/(1-?)))(y_{t}(j)-y_{t})

% 4. iterate to get future response
for i=2:T
    a(i)= rhoa*a(i-1)+ v(i);
    output(i)= rhoa*output(i-1) - sigma* psiya* (1-rhoa)* (1-beta*rhoa)*gamma*v(i);
    inflation(i)= rhoa*inflation(i-1) - sigma* psiya* (1-rhoa)* kappa *gamma*v(i);
    ir(i)=rho+ psipi*inflation(i) + psiy* output(i);
    y(i)=psiya*a(i)+ output(i);
    r(i)=ir(i)-inflation(i);
    price(i)= price(i-1)- output(i);
    product(i)= -epsilon*(price(i))+ output(i);
    mc(i)=  (sigma+(phi+alpha)/(1-alpha))*product(i);
end

% 3. Consider the level too
outputn= psiya.*a+ (1-alpha)*(0-log(1-alpha))/((1-alpha)*sigma+ phi+ alpha);
inflationn= zeros(T);
irn= rho-sigma*(1-rhoa)*psiya.*a;

% 4. Plot response to a one-time shock
t=100;
x =2:t;
subplot(3,3,1);
plot(x,a(2:t))
title('technology')
subplot(3,3,2);
plot(x,output(2:t))
title('output gap')
subplot(3,3,3);
plot(x,inflation(2:t))
title('inflation')
subplot(3,3,4);
plot(x,ir(2:t))
title('nominal rate')
subplot(3,3,5);
plot(x,r(2:t))
title('real rate')
subplot(3,3,6);
plot(x,y(2:t))
title('output')
subplot(3,3,7);
plot(x,price(2:t))
title('price deviation')
subplot(3,3,8);
plot(x,product(2:t))
title('product(j)')
subplot(3,3,9);
plot(x,mc(2:t))
title('marginal cost(j)')
suptitle('Economic path for random technology shocks')
print -djpeg -r600 q4_2.jpg

end
%% Change psipi and psiy in Taylor rule:
psipi= 10; % taylor rule of inflation
psiy= 0; % taylor rule of output

    
    
    
    
    