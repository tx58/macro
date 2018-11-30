%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Macro HW 9
% Tianli Xia November 28th
% This is the draft, use do instead.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
sigma= 2; % consumption elasiticity
phi= 3; % labor elasticity
epsilon= 5; % different goods bundle elasticity
alpha= 0.3; % production coefficient
beta= 0.99; % intertemporal discount factor
theta= 0.75; % ratio of firms which cannot change prices
psipi= 1.5; % taylor rule of inflation
psiy= 0.5; % taylor rule of output
rhoa= 0.8; % technology AR(1) coefficient
sigmaa2= 1; % technology AR(1) variance
rhoz= 0.5; % demand shock AR(1) coeffeicient
sigmaz2= 1; % demand shock AR(1) variance
rhov= 0; % policy shock AR(1) coefficient
sigmav2= 0.5; % policy shock AR(1) variance
rhou= 0.8; % cost-push shock AR(1) coefficient
sigmau2= 1; % cost-push shock AR(1) coefficient
period= 100;

rho= -log(beta);  % linearied beta
psiya= (1+phi)/(sigma*(1-alpha)+phi+alpha); % tech effect on output, affecting efficient level of output 
lamda= (1-theta)*(1-beta*theta)/theta*(1-alpha)/(1-alpha+alpha*epsilon);
kappa= lamda*(sigma+(phi+alpha)/(1-alpha));
gamma= 1/((1-beta*rhoa)*(sigma*(1-rhoa)+psiy)+ kappa*(psipi-rhoa));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stable/unstable Uncoulping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A0=[1,0,0,0,0,0;
    0,1,0,0,0,0;
    0,0,1,0,0,0;
    0,0,0,1,0,0;
    0,0,0,0,beta,0;
    0,0,0,0,1/sigma,1];
A1=[rhoa,0,0,0,0,0;
    0,rhoz,0,0,0,0;
    0,0,rhov,0,0,0;
    0,0,0,rhou,0,0;
    0,0,0,1,1,-kappa;
    (1-rhoa+psiy/sigma),-(1-rhoz)/sigma,1/sigma,0,psipi/sigma,1+psiy/sigma];
C1=[sqrt(diag([sigmaa2, sigmaz2, sigmav2, sigmau2])); zeros(2,4)];
A=A0\A1;
C=A0\C1;
C=C(1:4,1:4);

cutoff =1; % cutoff of stable/unstable variables
egen= abs(eig(A)) < cutoff;
n1=4; % number of predetermined variables; 
n2=2; % number of jump variables;
n= n1+n2; % total number of variables;

% Complex generalized Schur decomposition

[ZZ,TT]=schur(A); % Schur Decompostion to find A=UTU*
[Z,T] = ordschur(ZZ,TT,'UDI'); % Reorder to put unstable part on the RLS 

% Alternatively, as in the slide:
% [AA,BB,Q,Z] = qz(A,B) for square matrices A and B, produces upper quasitriangular matrices AA and BB, 
% and unitary matrices Q and Z such that Q*A*Z = AA, and Q*B*Z = BB. 
% For complex matrices, AA and BB are triangular.
% [S,T,Q,Z]= qz(eye(size(A)),A); % Q*A*Z= T, Q*S*Z= I. A=inv(Q)*T*inv(Z);
% [S,T,Q,Z]= reorder([S,T,Q,Z]);
% Reorder step would be problematic because it involves complex numbers

logcon = abs(diag(T)) <= cutoff ;

if sum(logcon) <n1
    warning("Too few stable roots: no stable solutions.");
    M= NaN; G=NaN; J0=NaN;
    return;
elseif sum(logcon) >n1
    warning("Too many stable roots: infinite number of stable solutions.");
    M= NaN; G=NaN; J0=NaN;
    return;
end

Z11=Z(1:n1,1:n1);
Z21=Z(n1+1:n,1:n1); % Divide Z into a block matrix [Z11; Z21]

M= real(A1(1:n1, 1:n1)); % X1(t+1)= MX1(t)+ e(t+1)
G= real(Z21/Z11); % X2(t+1)= GX1(t) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulating 200 periods
% X1(t+1)= MX1(t)+ Ce(t+1),  X1(t)=[a(t) z(t) v(t) u(t)]
% X2(t)= GX1(t) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step1. Impulse response
% 1. Specify the shock
shock= zeros(n1,period);
shock(1,2)= 1;
X1(:,1)= zeros(4,1); % initial shock level
a(1)=0;
z(1)=0;
v(1)=0;
u(1)=0;
% 2. Calculate the deviation of output and inflation
X2(:,1)= G*X1(:,1);
outputgap(1)= 0; % start from steady state
inflation(1)= 0;
ir(1)=rho;
r(1)= ir(1);

% 4. iterate to get future response
for i=2:period
    X1(:,i)= M*X1(:,i-1)+ C*shock(:,i);
    a(i)= X1(1,i);
    z(i)= X1(2,i);
    v(i)= X1(3,i);
    u(i)= X1(4,i);
    X2(:,i)= G*X1(:,i); 
    inflation(i)= X2(1,i); % pi
    outputgap(i)= X2(2,i); % x 
    ir(i)=rho+ psipi*inflation(i) + psiy* outputgap(i)+ v(i);
    y(i)=psiya*a(i)+ outputgap(i); % yt= yte+ x
    outputn(i)= y(i)- u(i)/kappa; % ytn=yte+ut/kappa
    r(i)=ir(i)-inflation(i); % r= i- pi
end

% 3. Consider the level too:
outpute= psiya.*a+ (1-alpha)*(0-log(1-alpha))/((1-alpha)*sigma+ phi+ alpha); % x=yt-yte -> yt=yte+x
irn= rho-sigma*(1-rhoa)*psiya.*a;

% 4. Plot response to a one-time shock
t=20;
x =2:t;
subplot(3,3,1);
plot(x,a(2:t))
title('technology')
subplot(3,3,2);
plot(x,z(2:t))
title('demand shock')
subplot(3,3,3);
plot(x,v(2:t))
title('monetary shock')
subplot(3,3,4);
plot(x,u(2:t))
title('cost shock')
subplot(3,3,5);
plot(x,outputgap(2:t))
title('output gap')
subplot(3,3,6);
plot(x,inflation(2:t))
title('inflation')
subplot(3,3,7);
plot(x,ir(2:t))
title('nominal rate')
subplot(3,3,8);
plot(x,r(2:t))
title('real rate')
subplot(3,3,9);
plot(x,y(2:t))
title('output')

suptitle('Impulse Response for a positive technology shock')
print -djpeg -r600 q4_1.jpg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kalmin Filter to recover 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step0. Speicfy the model
% Xtt(i)=M*Xtt(i-1)+ C*e(i-1)
% Ztt(i)=G*Xtt(i) +W*w(i)
W=zeros(2);
ztt=X2;
% Step1. recover latent variables
ptt1(:,:,1)= inv(eye(n1)-M)*C*C';
xtt1(:,1)= zeros(4,1);
i=1;
kt(:,:,i)=ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W');
ptt(:,:,i)=ptt1(:,:,i)-ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W')*G*ptt1(:,:,i);
ptt1(:,:,i+1)=M*ptt(:,:,i)*M'+C*C';
xtt1(:,i+1)= M*xtt1(:,i);
xtt(:,i)= M*zeros(4,1) +kt(:,:,i)*(ztt(:,i)-G*xtt1(:,i));

for i=2:period
    kt(:,:,i)=ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W');
    ptt(:,:,i)=ptt1(:,:,i)-ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W')*G*ptt1(:,:,i);
    ptt1(:,:,i+1)=M*ptt(:,:,i)*M'+C*C';
    xtt1(:,i+1)= M*xtt1(:,i);
    xtt(:,i)= M*xtt(:,i-1)+kt(:,:,i)*(ztt(:,i)-G*xtt1(:,i));
end

% Step2. Compare the estimated value and realized value
time=1:1:period;
subplot(2,2,1);
plot(time, xtt(1,:))
hold on
plot(time, a)
legend("Prediceted value","Realized value")
xlabel("Time")
title("Technology Shock")

subplot(2,2,2);
plot(time, xtt(2,:))
hold on
plot(time, z)
legend("Prediceted value","Realized value")
xlabel("Time")
title("Demand Shock")

subplot(2,2,3);
plot(time, xtt(3,:))
hold on
plot(time, v)
legend("Prediceted value","Realized value")
xlabel("Time")
title("Monetary Policy Shock")

subplot(2,2,4);
plot(time, xtt(4,:))
hold on
plot(time, u)
legend("Prediceted value","Realized value")
xlabel("Time")
title("Cost-push Shock")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alternative way to run Kalmin Filter (MATLAB command)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plant = ss(M,C,G,0,-1,'inputname','w','outputname','X2');
Q = eye(n1); 
R = zeros(2);
[kalmf,A,B,C] = kalman(Plant,Q,R);