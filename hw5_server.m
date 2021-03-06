% Tianli Xia
% Macro Homework 5
% Dynamic programming: value-iteration
% Value Function:
%% 
% $V(k, y)= max_{k'} log(e^{y_{t}}k^{\alpha} - (1-\delta)k - k') + \beta EV(k',y')$
% EV(k',y')=\Pi(y_{t+1}|y_{t})V(k_{t+1},y_{y+1})$
% Policy Function:
%% 
% $Q(k,y; k')= log(e^{y_{t}}k^{\alpha} - (1-\delta)k - k') + \beta EQ(k',y'; k'')$

%%
% Other important notations:
% 
% * Transition: k=k', P matrix
% * State variable: k, y
% * Action varible: k'
% 

clear all
%% Initialize the parameters

alpha=0.35; % the production rate: $f(k)=k^{\alpha}$
beta=0.95; % the intertemporary patience
delta=0.1; % capital deprecation rate
k=1000; % number of grids
len= 0.02; % length of grid, it is good to start from big value and then decrease.
% In pratice, I set length to be 0.01 at the beginning, at find that the
% optimal value at steady state should be within (0, 0.2), hence I set the
% length to a lower value to cover this interval.
start= len;
state= start:len:start+len*(k-1); % different states in grids:
action= start:len:start+len*(k-1);

%% Initialize the transiontion matrix P
m= 7; % number of discrete points we approximate
lamda= 0.98; % coefficient of AR(1) process
sigmaY= sqrt(0.1); % standard deviation of Y  
sigmaE= sqrt(1-lamda^2)*sigmaY; % standard deviation of Y
Y(m+2)=inf;
Y(1)= -inf; % Set the boundary value
P(m,m)= 0; % Define (P(t,t-1))
for i=1:m
    Y(i+1)=(i-((m+1)/2))*sigmaY;
end
for i=1:m % Note that in the loop i=1-> Y(i)=-inf, so we need P(1,.)=Y(2)
    for j=1:m
        P(i,j)=normcdf(((Y(j+1)+Y(j+2))/2-lamda*Y(i+1))/sigmaE)-...
               normcdf(((Y(j+1)+Y(j))/2-lamda*Y(i+1))/sigmaE); 
    end
end
    
%% Initialize value function matrix
V = ones(k,m); % V: state* action matrix
PI = ones(k,m); % PI: best policy matrix
threshold= 10^(-4); % Tolerance level in the loop
epsilon=1; % random initial value of the gao between two loops
count=0; % Count how many times the loop runs

tic
while epsilon> threshold % loop until $\epsilon$ converges
    V_temp=-inf*ones(k,m,k); % Any infeasible capital brings -inf value
    for s=1:k % value function iteration: get current value
        for j=1:m
            a_max=exp(Y(j+1))*state(s)^alpha+ (1-delta)*state(s); % Calculate feasible action sets: $0<=a<=state(s)^alpha+ (1-delta)*state(s)$
            for a=1:min(k, ceil( (a_max-start)/len) )
                V_temp(s,j,a)= log( exp(Y(j+1))*state(s)^alpha+ (1-delta)*state(s) -action(a)) + beta*P(j,:)*V(a,:)'; 
                % This directly comes from Bellman Optimality equation
            end
        end
        [V_new, PI]= max(V_temp, [], 3); % Calculate (1) New value function; (2) best action function.
    end
    epsilon= ( max(max( abs(V_new -V)))) % Calculate the current error
    V=V_new;
    count=count+1; % Count times the loop ends
end
toc 
%fprintf("The loop ends in %d runs, the gap within final two loops are %f.", count, epsilon)

save('solution.mat', 'V','PI', 'count', 'epsilon')
save('solution_1000.mat')