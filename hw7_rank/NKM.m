classdef NKM
   properties
      sigma= 2; % consumption elasiticity 
      phi= 3; % labor elasticity
      epsilon= 5; % different goods bundle elasticity
      alpha= 0.3; % production coefficient
      rhoa= 0.8; % technology AR(1) coefficient
      sigmaa2= 1; % technology AR(1) residuals
      beta= 0.99; % intertemporal discount factor
      theta= 0.75; % ratio of firms which cannot change prices
      psipi= 1.5; % taylor rule of inflation
      psiy= 0.5; % taylor rule of output
      T % total periods
      % Some definitions:
      rho % linearied beta
      psiya % tech effect on output
      lamda
      kappa
      gamma  
      
   end
   methods
       function obj = init(obj)
          obj.rho= -log(obj.beta);  % linearied beta
          obj.psiya= (1+obj.phi)/(obj.sigma*(1-obj.alpha)+obj.phi+obj.alpha); % tech effect on output
          obj.lamda= (1-obj.theta)*(1-obj.beta*obj.theta)/obj.theta*(1-obj.alpha)/(1-obj.alpha+obj.alpha*obj.epsilon);
          obj.kappa= obj.lamda*(obj.sigma+(obj.phi+obj.alpha)/(1-obj.alpha));
          obj.gamma= 1/((1-obj.beta*obj.rhoa)*(obj.sigma*(1-obj.rhoa)+obj.psiy)+ obj.kappa*(obj.psipi-obj.rhoa));
       end
       function r = shock(obj, v ,str, name)
           % 1. Specify the shock
           obj.T= size(v,1); 
           a(1)= 0; % initial technology level
           
           % 2. Calculate the deviation of output and inflation
           output(1)= 0; % start from steady state
           inflation(1)=0;
           ir(1)=obj.rho;
           r(1)= ir(1);
           
           % 3. price-related variables:
           price(1)= obj.theta*inflation(1)/(1-obj.theta); % price deviation from optimal p1*-pt
           product(1)= -obj.epsilon*(price(1))+ output(1); % output for product j: y_{t}(j)=-?(p_{t}^{j}-p_{t+k})+y_{t+k}
           mc(1)= (obj.sigma+(obj.phi+obj.alpha)/(1-obj.alpha))*product(1); % mc for product j: mc_{t}(j)-mc_{t}= (?+((?+?)/(1-?)))(y_{t}(j)-y_{t})
           
           % 4. iterate to get future response
           for i=2:obj.T
               a(i)= obj.rhoa*a(i-1)+ v(i);
               output(i)= obj.rhoa*output(i-1) - obj.sigma* obj.psiya* (1-obj.rhoa)* (1-obj.beta*obj.rhoa)*obj.gamma*v(i);
               inflation(i)= obj.rhoa*inflation(i-1) - obj.sigma* obj.psiya* (1-obj.rhoa)* obj.kappa *obj.gamma*v(i);
               ir(i)=obj.rho+ obj.psipi*inflation(i) + obj.psiy* output(i);
               y(i)=obj.psiya*a(i)+ output(i);
               r(i)=ir(i)-inflation(i);
               price(i)= price(i-1)- output(i);
               product(i)= -obj.epsilon*(price(i))+ output(i);
               mc(i)=  (obj.sigma+(obj.phi+obj.alpha)/(1-obj.alpha))*product(i);
           end
           
           % 3. Consider the level too
           outputn= obj.psiya.*a+ (1-obj.alpha)*(0-log(1-obj.alpha))/((1-obj.alpha)*obj.sigma+ obj.phi+ obj.alpha);
           inflationn= zeros(obj.T);
           irn= obj.rho-obj.sigma*(1-obj.rhoa)*obj.psiya.*a;
           
           % 4. Plot response to a one-time shock
           t=obj.T;
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
           suptitle(name)
           print(str,'-djpeg','-r600')
           snapnow
       end
   end
end


