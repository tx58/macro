classdef gd
   properties
       alpha=0.33; % the production rate: $f(k)=k^{\alpha}$
       beta=0.98; % the intertemporary patience
       delta=0.05; % capital deprecation rate
       N=600; % number of grids
       T=200; % number of periods
       len= 0.02; % length of grid, it is good to start from big value and then decrease.
       sigma= 0.8; % u(c)= c^(1-sigma)/(1-sigma)
       A= 1;
       % In pratice, I set length to be 0.01 at the beginning, at find that the
       % optimal value at steady state should be within (0, 0.2), hence I set the
       % length to a lower value to cover this interval.
       %start= len;
       %state= start:len:start+len*(N-1); % different states in grids:
       %action= start:len:start+len*(N-1);
       %terminal= start+len*(N-1);
       
       
       % Some definitions:
       k_star
       c_star
       lamda_star
       x_star
       c
       k
       x
       lamda
       k0
       V
       PI
       
       % Some results we get from the paras:
    
   end
   methods
       function obj = solvess(obj)
           obj.k_star= ((1/obj.beta- 1+ obj.delta )/ (obj.A*obj.alpha) )^(1/(obj.alpha-1) );
           obj.c_star= obj.A*obj.k_star^obj.alpha - obj.delta* obj.k_star;
           obj.lamda_star= obj.c_star ^ (-obj.sigma);
       end
       
       function obj = setstart(obj, k)
           obj.k0= k;
       end
       
       function obj = shooting(obj, tolerance)
           if nargin<2
               tolerance=1e-5;
           end
            obj.c=zeros(obj.T,1);
           obj.k=zeros(obj.T,1);
           obj.lamda=zeros(obj.T,1);
           obj.k(1)=obj.k0;
           cmax= obj.k(1)^obj.alpha+ obj.k(1)*(1-obj.delta);
           cmin= 0;
           count=1;
           while ( abs((obj.k(obj.T)-obj.k_star)/obj.k_star)>tolerance )
               c0 = 1/2*(cmax+cmin);
               obj.c(1)= c0;
               for t=2:obj.T
                   obj.k(t)= max(0, (1-obj.delta)*obj.k(t-1)+ ...
                       obj.A*obj.k(t-1)^obj.alpha - obj.c(t-1)) ;
                   obj.c(t)= max(0, obj.c(t-1)*(obj.beta*( (1-obj.delta) + ...
                       obj.A*obj.alpha*obj.k(t)^(obj.alpha-1)))^(1/obj.sigma)) ;
                   
               end
               if obj.k(obj.T) > obj.k_star && obj.c(obj.T) < obj.c_star
                   cmin= c0;
               else
                   cmax= c0;
               end
               abs(obj.k(obj.T)-obj.k_star)
               count=count+1;
           end
           fprintf("The loop ends in %d runs, the gap within final two loops are %f.",count, abs(obj.k(obj.T)-obj.k_star) )
           for t=1:199
               obj.x(t)= obj.A*obj.k(t)^obj.alpha- obj.c(t);
           end
           obj.x_star= obj.A*obj.k_star^obj.alpha- obj.c_star;
       end
       
       function obj = vf(obj, threshold)
           if nargin<2
               threshold=1e-5;% Tolerance level in the loop
           end

           obj.V = ones(obj.N,1); % V: state* action matrix
           obj.PI = ones(obj.N,1); % PI: best policy matrix
           epsilon=1; % random initial value of the gao between two loops
           count=0; % Count how many times the loop runs
           
           while epsilon> threshold % loop until $\epsilon$ converges
               V_temp=-inf*ones(obj.N,obj.N); % Any infeasible capital brings -inf value
               for s=1:obj.N % value function iteration: get current value
                   a_max= state(s)^obj.alpha+ (1-obj.delta)*state(s) ; % Calculate feasible action sets: $0<=a<=state(s)^alpha+ (1-delta)*state(s)$
                   for a=1:min(obj.N, ceil( (a_max-start)/obj.len) )
                       V_temp(s,a)= (state(s)^obj.alpha+ (1-obj.delta)*state(s) -action(a))^(1-obj.sigma)/(1-obj.sigma) ...
                           + obj.beta* obj.V(a);
                       % This directly comes from Bellman Optimality equation
                   end
                   [V_new, obj.PI]= max(V_temp, [], 2); % Calculate (1) New value function; (2) best action function.
               end
               epsilon= ( max( abs(V_new -obj.V))); % Calculate the current error
               obj.V=V_new;
               count=count+1; % Count times the loop ends
           end          
           fprintf("The loop ends in %d runs, the gap within final two loops are %f.", count, epsilon)       
       end
       
       function obj = vf_finite(obj, T, V2, threshold)
           if nargin<4
               threshold=1e-5;% Tolerance level in the loop
           end
           obj.V = ones(obj.N, T); % V: state* action matrix
           obj.PI = ones(obj.N, T); % PI: best policy matrix
           threshold= 10^(-4); % Tolerance level in the loop
           epsilon=1; % random initial value of the gap between two loops
           count=0; % Count how many times the loop runs

           while epsilon> threshold % loop until $\epsilon$ converges
               V_temp2=-inf*ones(obj.N,T,obj.N); % Any infeasible capital brings -inf value
               for s=1:obj.N % value function iteration: get current value
                   a_max=state(s)^obj.alpha+ (1-obj.delta)*state(s); % Calculate feasible action sets: $0<=a<=state(s)^alpha+ (1-delta)*state(s)$
                   for t=1:10
                       for a=1:min(obj.N, ceil( (a_max-start)/obj.len) )
                           if t~=10
                               V_temp2(s,t,a)= (state(s)^obj.alpha+ (1-obj.delta)*state(s) -action(a))^(1-obj.sigma)/(1-obj.sigma) ...
                                   + obj.beta* V2(a, t+1);
                               % This directly comes from Bellman Optimality equation
                           else
                               V_temp2(s,t,a)= (state(s)^obj.alpha+ (1-obj.delta)*state(s) -action(a))^(1-obj.sigma)/(1-obj.sigma) ...
                                   + obj.beta* obj.V(a);
                           end
                       end
                   end
                   [V_new2, obj.PI]= max(V_temp2, [], 3); % Calculate (1) New value function; (2) best action function.
               end
               epsilon= ( max(max( abs(V_new2 -V2)))); % Calculate the current error
               V2=V_new2;
               count=count+1; % Count times the loop ends
           end
           
           fprintf("The loop ends in %d runs, the gap within final two loops are %f.", count, epsilon)
    
       end
       
       function lamda= findpath(obj, k0)
           obj.c=zeros(obj.T,1);
           obj.k=zeros(obj.T,1);
           obj.lamda=zeros(obj.T,1);
           obj.k(1)=k0;
           cmax= obj.k(1)^obj.alpha+ obj.k(1)*(1-obj.delta);
           cmin= 0;
           count=1;
           while ( abs((obj.k(obj.T)-obj.k_star)/obj.k_star)>10e-4 )
               c0 = 1/2*(cmax+cmin);
               obj.c(1)= c0;
               for t=2:obj.T
                   obj.k(t)= max(0, (1-obj.delta)*obj.k(t-1)+ ...
                       obj.A*obj.k(t-1)^obj.alpha - obj.c(t-1)) ;
                   obj.c(t)= max(0, obj.c(t-1)*(obj.beta*( (1-obj.delta) + ...
                       obj.A*obj.alpha*obj.k(t)^(obj.alpha-1)))^(1/obj.sigma)) ;
                   
               end
               if obj.k(obj.T) > obj.k_star && obj.c(obj.T) < obj.c_star
                   cmin= c0;
               else
                   cmax= c0;
               end
               abs(obj.k(obj.T)-obj.k_star)
               count=count+1;
           end
           lamda= c0^(-obj.sigma);
       end
       
   end
  
   
end


