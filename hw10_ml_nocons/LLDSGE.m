function [LL]=LLDSGE(THETA)
global Z

tau       = THETA(1);  % Since I am using sigma for the standard deviaion of the shocks, I am using tau to denote the CRRA parameter.
beta      = THETA(2); % discount factor
theta     = THETA(3); % degree of price stickiness
phi_pi    = THETA(4); % taylor rule parameter
phi_y     = THETA(5); % taylor rule parameter
varphi    = THETA(6); %inverse of elastiicity of labor supply
alpha     = THETA(7); %production function parameter
eps       = THETA(8); % elasticity of substitution between goods i and j in the consumption basket
rho_v     = THETA(9); %persistence parameter
rho_a     = THETA(10); %persistence parameter
rho_z     = THETA(11); %persistence parameter
rho_u     = THETA(12); %persistence parameter
sigma_v   = THETA(13); %standard deviation
sigma_a   = THETA(14); %standard deviation of innovation to a_t
sigma_z   = THETA(15); %srandard deviation of z
sigma_u   = THETA(16); %standard deviation of u

T = DSGE_solve(tau,beta,theta,phi_pi,phi_y,varphi,alpha,eps,rho_v,rho_a,rho_z,rho_u,sigma_v,sigma_a,sigma_z,sigma_u);
% Return the solved State space function parameters
psi_yna = (1+varphi)/(tau*(1-alpha)+varphi+alpha);
ye= (1-alpha)*(0-log(1-alpha))/((1-alpha)*tau+ varphi+ alpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put model in state space form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=diag([rho_a,rho_z,rho_v,rho_u]);
C=diag([sigma_a,sigma_z,sigma_v,sigma_u]);
D=[T; zeros(2,4)] + [zeros(1,4); psi_yna 0 0 0; 0 0 1 0; -1/(1-alpha) 0 0 0 ];
temp=Z;
% Realign measurement equation
ZZ(1,:)=temp(2,:); % inflation
ZZ(2,:)=temp(1,:); % output
ZZ(3,:)=temp(3,:)- phi_pi*temp(2,:) - phi_y*temp(1,:); % interest rate
ZZ(4,:)=temp(4,:)- temp(1,:)/(1-alpha); % labor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for Kalman filter etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
period=length(Z);
LL=0;
P= dlyap(A,C*C');%Initial uncertainty equal to unconditional variance of state
Xfilt=zeros(4,1); %initial value for the filtered state.
CC=C*C';
%RR=zeros(4,4);
RR=eye(4)*.001;

%Compute recursive likelihood using the Kalman filter
for tt=1:period
    a=ZZ(1:4,tt)-D*Xfilt;%These are the innovations (i.e. Ztilde)
    Omega=D*P*D'+RR;
    Omegainv=eye(4)/(Omega);
    K=P*D'*Omegainv;
    Xfilt=A*Xfilt+A*K*a;
    P = A*(P-P*D'*Omegainv*D*P)*A' + CC;
    LL = LL - 0.5*(log(det(Omega)) + a'*Omegainv*a);
end
% if min(eu)==0
%     LL=-9e+200;
% end
% if imag(LL)~=0
%     LL=-9e+200;
% end

LL=-LL;%Because simulated annealing is minimizing
