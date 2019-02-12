% Tianli Xia
% Macro Homework 3
% Dynamic programming: value-iteration
% Policy Function:
%% 
% $Q(s, k)= log(s^{\alpha} - (1-\delta)s - k) + \beta V(k)$

% Value Function:
%% 
% $V(s)= max_{k} Q(s, k)$

%%
% Other important notations:
% 
% * Transition: s=k
% * State variable: 
% * Action varible: k
% 

clear all

%% Initialize the parameters

A=1; % techonology
alpha=0.33; % the production rate: $f(k)=k^{\alpha}$
beta=0.98; % the intertemporary patience
delta=0.05; % capital deprecation rate
k=600; % number of grids
len= 0.002; % length of grid, it is good to start from big value and then decrease.
sigma= 0.8; % u(c)= c^(1-sigma)/(1-sigma)
% In pratice, I set length to be 0.01 at the beginning, at find that the
% optimal value at steady state should be within (0, 0.2), hence I set the
% length to a lower value to cover this interval.
start= 10;
state= start:len:start+len*(N-1); % different states in grids:
action= start:len:start+len*(N-1);
terminal= start+len*(N-1);


%% Show the capital, consumption along the path
notpatient=gd; 
notpatient= notpatient.solvess;
notpatient= notpatient.setstart(0.85*notpatient.k_star);
notpatient= notpatient.shooting(1e-5);

patient=notpatient;
patient.beta=0.99;
patient= patient.solvess;
patient= patient.setstart(0.85*patient.k_star);
patient= patient.shooting;

hightech=notpatient;
hightech.A=1.01;
hightech= hightech.solvess;
hightech= hightech.setstart(0.85*hightech.k_star);
hightech= hightech.shooting;

t=100;
subplot(3,2,1)
plot(1:1:t, notpatient.k(1:t), 'b-', 1:1:t, ones(t,1)*notpatient.k_star, 'r--' )
title("capital")
subplot(3,2,3)
plot(1:1:t, notpatient.c(1:t), 'b-', 1:1:t, ones(t,1)*notpatient.c_star, 'r--' )
title("consumption")
subplot(3,2,5)
plot(1:1:t, notpatient.x(1:t), 'b-', 1:1:t, ones(t,1)*notpatient.x_star, 'r--' )
title("investment")
xlabel("period")
subplot(3,2,2)
plot(1:1:t, patient.k(1:t), 'b-', 1:1:t, ones(t,1)*patient.k_star, 'r--' )
title("capital")
subplot(3,2,4)
plot(1:1:t, patient.c(1:t), 'b-', 1:1:t, ones(t,1)*patient.c_star, 'r--' )
title("consumption")
subplot(3,2,6)
plot(1:1:t, patient.x(1:t), 'b-', 1:1:t, ones(t,1)*patient.x_star, 'r--' )
title("investment")
xlabel("period")
    print -djpeg -r600 hw3_patient

%% Plot optimal path

for i=1:1:200
    ki= 8.5+i/50;
    lamda1(i)=(-delta*ki+A*ki^alpha)^-sigma;
    lamda2(i)= ((1-delta)*ki+A*ki^alpha-k_star)^-sigma;
end
    plot(1:1:200, lamda1)
    hold on
    plot(1:1:200, lamda2)
    label("Feasibility","Eular")
    title("phase diagram")





