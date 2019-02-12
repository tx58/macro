% Set up and estimate miniture DSGE model
clc
clear all
close all
global Z
load('Z'); %load data. Order of variables: Inflation, output, interest rate, labor supply

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial values of structural parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau       = 3;  % Since I am using sigma for the standard deviaion of the shocks, I am using tau to denote the CRRA parameter.
beta      = 0.99; % discount factor
theta     = 3/4; % degree of price stickiness
phi_pi    = 1.5; % taylor rule parameter
phi_y     = 0.125; % taylor rule parameter
varphi    = 3; %inverse of elastiicity of labor supply
alpha     = 1/3; %production function parameter
eps       = 6; % elasticity of substitution between goods i and j in the consumption basket
rho_v     = 0.5; %persistence parameter
rho_a     = 0.75; %persistence parameter
rho_z     = 0.5; %persistence parameter
rho_u     = 0.3; %persistence parameter
sigma_v   = 0.02; %standard deviation
sigma_a   = (0.012^2*(1-rho_a^2))^.5; %standard deviation of innovation to a_t
sigma_z   = 0.01; %standard deviation
sigma_u   = 0.1; %standard deviation


THETA=[tau,beta,theta,phi_pi,phi_y,varphi,alpha,eps,rho_v,rho_a,rho_z,rho_u,sigma_v,sigma_a,sigma_z,sigma_u]';%Starting value for parameter vector
LB=[0,0,0,1,0,1,0,1,zeros(1,8)]';%Lower bound for parameter vector
UB=[10,1,1,5,5,10,1,50,1,1,1,1,10,10,10,10]';%Upper bound for parameter vector
x=THETA;

% x=THETA;
sa_t= 5;
sa_rt=.3;
sa_nt=5;
sa_ns=5;
% warning off all;

[xhat]=simannb( 'LLDSGE', x, LB, UB, sa_t, sa_rt, sa_nt, sa_ns, 1);

%--------------------------------------------------------------------------
thetalabel=['tau    ';'beta   ';'theta  ';'phi_pi ';'phi_y  ';'varphi ';'alpha  ';'eps    ';'rho_v  ';'rho_a  ';'rho_z  ';'rho_u  ';'sigmav ';'sigmaa ';'sigmaz ';'sigmau '];
disp('ML estimate of THETA')
disp([thetalabel, num2str(xhat)])
%--------------------------------------------------------------------------
% a_hat = filter_gap_DSGE(xhat,Z);
% figure(2)
% plot(Z(1,:),'linewidth',2);hold on;plot(a_hat,'linewidth',2)
% legend('Output','Output gap')

%% plot the impulse response function and predicted value
Model= DSGE;
Model= Model.setvalue(xhat);
Model= Model.init;
Model= Model.solve;
v=zeros(4,20);
Model.shock(v, 'q1_1', 'Steady state without positive technology shocks');
v(1,2)=0.1; % Set a technology shock 
Model.shock(v, 'q1_2', 'Impulse Response for a positive technology shock');
v=zeros(4,20);
v(2,2)=0.1; % Set a demand shock 
Model.shock(v, 'q1_3', 'Impulse Response for a positive demand shock');
v=zeros(4,20);
v(3,2)=0.1; % Set a monerary policy shock 
Model.shock(v, 'q1_4', 'Impulse Response for a positive monerary shock');
v=zeros(4,20);
v(4,2)=0.1; % Set a cost-push shock 
Model.shock(v, 'q1_5', 'Impulse Response for a positive cost-push shock');

%% Plot predicted value under estimated parameter
pi= Z(2,:);
y= Z(1,:);
ir= Z(3,:);
labor= Z(4,:);
%ytn= Model.psiya*a + Model.ye - u/Model.kappa;
ztt= [pi; y; ir-Model.psipi*pi-Model.psiy*(y); labor-y/(1-alpha)];
% Above is some transformation to make filter feasible
C= Model.C;
M= Model.M;
G= [Model.G; zeros(2,4)]+ [0 0 0 0; Model.psiya 0 0 0; 0 0 1 0; -1/(1-alpha) 0 0 0];
n1=4;

W=eye(4)*0.001;

% Step1. recover latent variables
ptt1(:,:,1)= dlyap(M,C*C');
xtt1(:,1)= zeros(4,1);
i=1;
kt(:,:,i)=ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W');
ptt(:,:,i)=ptt1(:,:,i)-ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W')*G*ptt1(:,:,i);
ptt1(:,:,i+1)=M*ptt(:,:,i)*M'+C*C';
xtt1(:,i+1)= M*xtt1(:,i);
xtt(:,i)= M*zeros(4,1) +kt(:,:,i)*(ztt(:,i)-G*xtt1(:,i));

Model.period=length(Z);
for i=2:Model.period
    kt(:,:,i)=ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W');
    ptt(:,:,i)=ptt1(:,:,i)-ptt1(:,:,i)*G'*inv(G*ptt1(:,:,i)*G'+W*W')*G*ptt1(:,:,i);
    ptt1(:,:,i+1)=M*ptt(:,:,i)*M'+C*C';
    xtt1(:,i+1)= M*xtt1(:,i);
    xtt(:,i)= M*xtt(:,i-1)+kt(:,:,i)*(ztt(:,i)-G*xtt1(:,i));
end

% Step2. use latent variable to find ytn
zthat= G*xtt;
% pihat= zthat(1,:);
% yhat=zthat(2,:)+ ye;
% irhat=zthat(3,:)- (-Model.psipi*pi-Model.psiy*(y-Model.ye)-Model.rho) ;
% laborhat=zthat(4,:) + y/(1-alpha);

pihat= zthat(1,:);
yhat=zthat(2,:);
irhat=zthat(3,:)+ Model.psipi*pi+ Model.psiy*(y);
laborhat=zthat(4,:) + labor-y/(1-alpha);


% Step3. Plot them
subplot(2,2,1)
time=1:1:Model.period;
plot(time, y)
hold on
plot(time, yhat)
legend("Realized value","Prediceted value")
xlabel("Time")
title("Output")
subplot(2,2,2)
time=1:1:Model.period;
plot(time, pi)
hold on
plot(time, pihat)
legend("Realized value","Prediceted value")
xlabel("Time")
title("Inflation")
subplot(2,2,3)
time=1:1:Model.period;
plot(time, ir)
hold on
plot(time, irhat)
legend("Realized value","Prediceted value")
xlabel("Time")
title("Interest rate")
subplot(2,2,4)
time=1:1:Model.period;
plot(time, labor)
hold on
plot(time, laborhat)
legend("Realized value","Prediceted value")
xlabel("Time")
title("Employment")
print -djpeg -r600 q2.jpg


