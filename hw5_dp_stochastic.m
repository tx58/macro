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
k=100; % number of grids
len= 0.2; % length of grid, it is good to start from big value and then decrease.
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
fprintf("The loop ends in %d runs, the gap within final two loops are %f.", count, epsilon)
save("solution.mat", "V","PI")

%% Plot graph
%% Plot the value function
% $V(k,y)= \max_{k'} log(exp^{y}s^{\alpha} - (1-\delta)s - k') + \beta P(y'|y) V(k',y')$
if exist('solution_1000.mat', 'file')
    load solution_1000.mat
end

for i=1:7
    plot(state, V(:,i))
    hold on
end
    grid on
    axis on
    xlabel("Current capital (k)")
    ylabel("Value (V(k,y))")
    lgd= legend("-0.9487","-0.6325","-0.3162","0","0.3162","0.6325","0.9487",'Location','southeast');
    title(lgd, "Current State (y)")
    title("Value Function on each state")
    hold off
                            
% print -djpeg -r600 hw5_value
%% Plot action function
%This plots show the optimal capital decision. At any point on the left side
%of steady state, we increase the capital and vice versa. At steady state
%it is stable. Hence the steady state is the crosspoint of the policy
%function and the 45 degree line.


for j=1:m
    i(j)=1;
    while i(j)~=PI(i(j),j) % By definition, this is the optimal decision
        i(j)=i(j)+1;
    end

    %fprintf("The optimal capital at the steady state is %f\n. The value at optimal capital is %f\n" ...
    %    , i*len, V(i,j) )
    %fprintf("If feasible, the planner will always choose steady state capital in the next period.")
        plot(state, action(PI(:,j)))
        x=state;
        y=state;
        hold on
        plot(x(i(j)),y(i(j)),'r*')
        
    %   print -djpeg -r600 hw5_action
end
    plot(x,y)
    hold on
    grid on
    axis on
    xlabel("Current capital")
    ylabel("Next period best action")
    title("Policy Function on each state")
    axis square
    hold off
    
    %% Plot potential function of actions in terms of action at each state 
% $Q(s,y,a)= \max_{k} log(exp^{y}s^{\alpha} - (1-\delta)s - k) + \beta V(k)$
for i=11:20:91
    y= reshape(V_temp(i,4,:), size(state));
    plot(state, y)
    hold on
end
xlabel("Next period capital decision")
ylabel("Optimality function q(s,a)")
lgd= legend("0.042","0.082","0.121","0.161","0.201","0.241",'Location','southeast');
title(lgd, "Current State (Captial)")
axis on
grid on
hold off

%% Calcualte the model predicted output, capital, and investment
% Start from a steady state, suppose all the decision is optimal (according to PI matrix)
% Simulate the process for 100 times

steady=repmat(1:1:k, m, 1); 
index= ceil ((sum((steady.*(steady==PI')), 2)./ sum(((steady==PI')),2 )));  % By definition, this is the optimal decision
% Here index denotes the steady state, we set it as the original value

r=100;
y(1)=4; % state variable y
start=index(y(1));
capital(150)=0; 
output(150)=0;
investment(150)=0;
capital(1)=state(start); % state variable capital
output(1)= exp(Y(y(1)+1))*capital(1)^alpha; % output according to production function
investment(1)=  action(PI(start,y(1)))- (1-delta)*capital(1); % investment decision 


for i=2:150
    y(i)=sum( rand(1)>=( cumsum(P(y(i-1),:))))+1; % random process of state 
    capital(i)= action(PI(index(y(i-1)),y(i-1))); % capital accumulation process
    output(i)= exp(Y(y(i)+1))*capital(i)^alpha;
    investment(i)=  action(PI(index(y(i)),y(i)))- (1-delta)*capital(i);
end

plot(1:i, investment)
hold on
plot(1:i, output)
hold on
plot(1:i, y)
