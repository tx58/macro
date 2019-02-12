clear all

%% Initialize the parameters

A=1; % techonology
alpha=0.33; % the production rate: $f(k)=k^{\alpha}$
beta=0.98; % the intertemporary patience
delta=0.05; % capital deprecation rate
N=600; % number of grids
k=N;
len= 0.0003; % length of grid, it is good to start from big value and then decrease.
sigma= 0.8; % u(c)= c^(1-sigma)/(1-sigma)
% In pratice, I set length to be 0.01 at the beginning, at find that the
% optimal value at steady state should be within (0, 0.2), hence I set the
% length to a lower value to cover this interval.
start= 10.01;
state= start:len:start+len*(N-1); % different states in grids:
action= start:len:start+len*(N-1); 
terminal= start+len*(N-1); 


%% Initialize value function matrix

A=1;
V = ones(k,1); % V: state* action matrix
PI = ones(k,1); % PI: best policy matrix
threshold= 10^(-4); % Tolerance level in the loop
epsilon=1; % random initial value of the gao between two loops
count=0; % Count how many times the loop runs

while epsilon> threshold % loop until $\epsilon$ converges
    V_temp=-inf*ones(k,k); % Any infeasible capital brings -inf value
    for s=1:k % value function iteration: get current value
        a_max= A*state(s)^alpha+ (1-delta)*state(s) ; % Calculate feasible action sets: $0<=a<=state(s)^alpha+ (1-delta)*state(s)$
        for a=1:min(k, ceil( (a_max-start)/len) )
            V_temp(s,a)= (A*state(s)^alpha+ (1-delta)*state(s) -action(a))^(1-sigma)/(1-sigma) ...
                + beta* V(a); 
            % This directly comes from Bellman Optimality equation
        end
        [V_new, PI]= max(V_temp, [], 2); % Calculate (1) New value function; (2) best action function.
    end
    epsilon= ( max( abs(V_new -V))) % Calculate the current error
    V=V_new;
    count=count+1; % Count times the loop ends
end

fprintf("The loop ends in %d runs, the gap within final two loops are %f.", count, epsilon)

V0= V;
PI0= PI;


%% Plot the value function
% $V(k)= \max_{k} log(s^{\alpha} - (1-\delta)s - k) + \beta V(k)$
plot(state, V)
grid on
axis on
xlabel("Current capital")
ylabel("Value")
title("Value Function on each state")
 
print -djpeg -r600 hw3_value
 
%% Plot potential function of actions in terms of action at each state 
% $Q(s,a)= \max_{k} log(s^{\alpha} - (1-\delta)s - k) + \beta V(k)$
for i=1:30:151
    plot(state, V_temp(i,:))
    hold on
end
xlabel("Next period capital decision")
ylabel("Optimality function q(s,a)")
lgd= legend("0.001","0.031","0.061","0.091","0.121","0.151",'Location','southeast');
title(lgd, "Current State (Captial)")
axis on
grid on
hold off
% for i=1:5:26
%     plot(state(PI(i)), V_temp(i,PI(i)), "r.")
%     hold on
% end
% hold off
%print -djpeg -r600 -hw3_control_2


%% Plot action function
%This plots show the optimal capital decision. At any point on the left side
%of steady state, we increase the capital and vice versa. At steady state
%it is stable. Hence the steady state is the crosspoint of the policy
%function and the 45 degree line.
i=1;
while i~=PI(i) % By definition, this is the optimal decision
    i=i+1;
end
      
fprintf("The optimal capital at the steady state is %f\n. The value at optimal capital is %f\n" ...
    , i*len, V(i) )
fprintf("If feasible, the planner will always choose steady state capital in the next period.")
    plot(state, action(PI))
    x=state;
    y=state;
    hold on
    plot(x,y)
    hold on
    plot(x(i),y(i),'r*')
    grid on
    axis on
    xlabel("Current capital")
    ylabel("Next period best action")
    title("Policy Function on each state")
    axis square
    print -djpeg -r600 hw3_action_2

%% New value function: given higher technology A=1.01    

A=1.01;
V = ones(k,1); % V: state* action matrix
PI = ones(k,1); % PI: best policy matrix
threshold= 10^(-4); % Tolerance level in the loop
epsilon=1; % random initial value of the gao between two loops
count=0; % Count how many times the loop runs

while epsilon> threshold % loop until $\epsilon$ converges
    V_temp=-inf*ones(k,k); % Any infeasible capital brings -inf value
    for s=1:k % value function iteration: get current value
        a_max= A*state(s)^alpha+ (1-delta)*state(s) ; % Calculate feasible action sets: $0<=a<=state(s)^alpha+ (1-delta)*state(s)$
        for a=1:min(k, ceil( (a_max-start)/len) )
            V_temp(s,a)= (A*state(s)^alpha+ (1-delta)*state(s) -action(a))^(1-sigma)/(1-sigma) ...
                + beta* V(a); 
            % This directly comes from Bellman Optimality equation
        end
        [V_new, PI]= max(V_temp, [], 2); % Calculate (1) New value function; (2) best action function.
    end
    epsilon= ( max( abs(V_new -V))) % Calculate the current error
    V=V_new;
    count=count+1; % Count times the loop ends
end

fprintf("The loop ends in %d runs, the gap within final two loops are %f.", count, epsilon)


%% Case2: consider increase techonology in period 10:
A=1;
T=10; % number of peiords
V2 = ones(k, T); % V: state* action matrix
PI2 = ones(k, T); % PI: best policy matrix
threshold= 10^(-4); % Tolerance level in the loop
epsilon=1; % random initial value of the gap between two loops
count=0; % Count how many times the loop runs

while epsilon> threshold % loop until $\epsilon$ converges
    V_temp2=-inf*ones(k,T,k); % Any infeasible capital brings -inf value
    for s=1:k % value function iteration: get current value
        a_max=A*state(s)^alpha+ (1-delta)*state(s); % Calculate feasible action sets: $0<=a<=state(s)^alpha+ (1-delta)*state(s)$
        for t=1:10
            for a=1:min(k, ceil( (a_max-start)/len) )
                if t~=10
                    V_temp2(s,t,a)= (A*state(s)^alpha+ (1-delta)*state(s) -action(a))^(1-sigma)/(1-sigma) ...
                        + beta* V2(a, t+1); 
                    % This directly comes from Bellman Optimality equation
                else
                    V_temp2(s,t,a)= (A*state(s)^alpha+ (1-delta)*state(s) -action(a))^(1-sigma)/(1-sigma) ...
                    + beta* V(a);
                end
            end
        end
        [V_new2, PI2]= max(V_temp2, [], 3); % Calculate (1) New value function; (2) best action function.
    end
    epsilon= ( max(max( abs(V_new2 -V2)))) % Calculate the current error
    V2=V_new2;
    count=count+1; % Count times the loop ends
end

fprintf("The loop ends in %d runs, the gap within final two loops are %f.", count, epsilon)


%% Plot graph

k0=find(PI0'==1:1:N, 1 );
k_path1(1)= k0;
k_path2(1)= k0; % start from current steady state
for t=1:60
    k_path1(t+1) = PI(k_path1(t));
    if t<10
        k_path2(t+1) = PI2(k_path2(t),t);
    else
        k_path2(t+1) = PI(k_path2(t));
    end    
end

k_path1= (k_path1-1)*len+ start;
k_path2= (k_path2-1)*len+ start;
c_path1= (1-delta)*k_path1(t)+ 1.01*k_path1(t)^alpha - k_path1(t+1);
c_path2= (1-delta)*k_path1(t)+ 1.01*k_path1(t)^alpha - k_path1(t+1);
x_path1= k_path1(t+1)- c_path1;
x_path2= 1.01*k_path2(t+1)^alpha- c_path2;

for t=1:50    
    c_path1(t+1) = (1-delta)*k_path1(t)+ 1.01*k_path1(t)^alpha - k_path1(t+1);
    x_path1(t+1) =  1.01*k_path1(t+1)^alpha- c_path1(t+1);

    if t<10
        c_path2(t+1) = (1-delta)*k_path1(t)+ k_path1(t)^alpha - k_path1(t+1);
        x_path2(t+1) = k_path2(t+1)^alpha - c_path2(t+1);
    else
        c_path2(t+1) = (1-delta)*k_path1(t)+ 1.01*k_path1(t)^alpha - k_path1(t+1);
        x_path2(t+1) = 1.01*k_path2(t+1)^alpha- c_path2(t+1);
    end
end
subplot(3,2,1)
plot(2:1:t, k_path1(2:t), 'b-', 2:1:t, ones(t-1,1)*k_path1(t), 'r--' )
title("capital")
subplot(3,2,3)
plot(2:1:t, c_path1(2:t), 'b-', 2:1:t, ones(t-1,1)*c_path1(t), 'r--' )
title("consumption")
subplot(3,2,5)
plot(2:1:t, x_path1(2:t), 'b-', 2:1:t, ones(t-1,1)*x_path1(t), 'r--' )
title("investment")
xlabel("period")
subplot(3,2,2)
plot(2:1:t, k_path2(2:t), 'b-', 2:1:t, ones(t-1,1)*k_path1(t), 'r--' )
title("capital")
subplot(3,2,4)
plot(2:1:t, c_path2(2:t), 'b-', 2:1:t, ones(t-1,1)*c_path1(t), 'r--' )
title("consumption")
subplot(3,2,6)
plot(2:1:t, x_path2(2:t), 'b-', 2:1:t, ones(t-1,1)*x_path1(t), 'r--' )
title("investment")
xlabel("period")
print -djpeg -r600 hw3_transition

