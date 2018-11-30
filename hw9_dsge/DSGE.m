classdef DSGE
   properties
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
      period % total periods
      
      % Some definitions:
      rho % linearied beta
      psiya % tech effect on output
      lamda
      kappa
      gamma  
      
      % Some results we get from the paras:
      C % covaraince matrix of shocks
      M % transition coefficient of latent vars
      G % transition coefficient of obeserved vars
      X2 % Observed sequence of variables (pi, output_gap, [interest rate])
      ye % efficient level output
      yte
   end
   methods
       function obj = init(obj)
          obj.rho= -log(obj.beta);  % linearied beta
          obj.psiya= (1+obj.phi)/(obj.sigma*(1-obj.alpha)+obj.phi+obj.alpha); % tech effect on output
          obj.lamda= (1-obj.theta)*(1-obj.beta*obj.theta)/obj.theta*(1-obj.alpha)/(1-obj.alpha+obj.alpha*obj.epsilon);
          obj.kappa= obj.lamda*(obj.sigma+(obj.phi+obj.alpha)/(1-obj.alpha));
          obj.gamma= 1/((1-obj.beta*obj.rhoa)*(obj.sigma*(1-obj.rhoa)+obj.psiy)+ obj.kappa*(obj.psipi-obj.rhoa));
          obj.ye= (1-obj.alpha)*(0-log(1-obj.alpha))/...
               ((1-obj.alpha)*obj.sigma+ obj.phi+ obj.alpha);
       end
       function obj = solve(obj) % Numerically find the solution
          A0=[1,0,0,0,0,0;
              0,1,0,0,0,0;
              0,0,1,0,0,0;
              0,0,0,1,0,0;
              0,0,0,0,obj.beta,0;
              0,0,0,0,1/obj.sigma,1];
          A1=[obj.rhoa,0,0,0,0,0;
              0,obj.rhoz,0,0,0,0;
              0,0,obj.rhov,0,0,0;
              0,0,0,obj.rhou,0,0;
              0,0,0,1,1,-obj.kappa;
              (1-obj.rhoa+obj.psiy/obj.sigma),-(1-obj.rhoz)/obj.sigma,1/obj.sigma,0,obj.psipi/obj.sigma,1+obj.psiy/obj.sigma];
          C1=[sqrt(diag([obj.sigmaa2, obj.sigmaz2, obj.sigmav2, obj.sigmau2])); zeros(2,4)];
          A=A0\A1;
          obj.C=A0\C1;
          obj.C=obj.C(1:4,1:4);
          
          cutoff =1; % cutoff of stable/unstable variables
          egen= abs(eig(A)) < cutoff;
          n1=4; % number of predetermined variables;
          n2=2; % number of jump variables;
          n= n1+n2; % total number of variables;
          
          % Complex generalized Schur decomposition
          
          [ZZ,TT]=schur(A); % Schur Decompostion to find A=UTU*
          [Z,T] = ordschur(ZZ,TT,'UDI'); % Reorder to put unstable part on the RLS          
          logcon = abs(diag(T)) <= cutoff ;
          
          if sum(logcon) <n1
              warning("Too few stable roots: no stable solutions.");
              obj.M= NaN; obj.G=NaN; J0=NaN;
          elseif sum(logcon) >n1
              warning("Too many stable roots: infinite number of stable solutions.");
              obj.M= NaN; obj.G=NaN; J0=NaN;
          end
          
          Z11=Z(1:n1,1:n1);
          Z21=Z(n1+1:n,1:n1); % Divide Z into a block matrix [Z11; Z21]
          
          obj.M= real(A1(1:n1, 1:n1)); % X1(t+1)= MX1(t)+ e(t+1)
          obj.G= real(Z21/Z11); % X2(t+1)= GX1(t)
       end
       
       function [X1, X2, ir, y, r] = shock(obj, shock ,str, name)
           % 1. Specify the shock
           obj.period= size(shock,2); 
           
           X1(:,1)= zeros(4,1); % initial shock level
           a(1)=0;
           z(1)=0;
           v(1)=0;
           u(1)=0;
           % 2. Calculate the deviation of output and inflation
           X2(:,1)= obj.G*X1(:,1);
           outputgap(1)= 0; % start from steady state
           inflation(1)= 0;
           ir(1)=obj.rho;
           r(1)= ir(1);
           
           % 4. iterate to get future response
           for i=2:obj.period
               X1(:,i)= obj.M*X1(:,i-1)+ obj.C*shock(:,i);
               a(i)= X1(1,i);
               z(i)= X1(2,i);
               v(i)= X1(3,i);
               u(i)= X1(4,i);
               X2(:,i)= obj.G*X1(:,i);
               inflation(i)= X2(1,i); % pi
               outputgap(i)= X2(2,i); % x
               ir(i)=obj.rho+ obj.psipi*inflation(i) + obj.psiy* outputgap(i)+ v(i);
               y(i)=obj.psiya*a(i)+ outputgap(i); % yt= yte+ x
               outputn(i)= y(i)- u(i)/obj.kappa; % ytn=yte+ut/kappa
               r(i)=ir(i)-inflation(i); % r= i- pi
           end
           
           % 3. Consider the level too:
           obj.yte= obj.psiya.*a+ (1-obj.alpha)*(0-log(1-obj.alpha))/...
               ((1-obj.alpha)*obj.sigma+ obj.phi+ obj.alpha); % x=yt-yte -> yt=yte+x
           irn= obj.rho-obj.sigma*(1-obj.rhoa)*obj.psiya.*a;
           
           % 4. Plot response to the shock
           t=obj.period;
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
           suptitle(name)
           print(str,'-djpeg','-r600')
           snapnow
       end
       
       function output = kalmin(obj, X1, X2,  str, name) % Use Kalmin filter to form estimates for latent vars
           if nargin > 3
               
           else
               name="Kalmin filter on recovering shocks";
               str="q3";
           end
           % Step0. Speicfy the model
           % Xtt(i)=M*Xtt(i-1)+ C*e(i-1)
           % Ztt(i)=G*Xtt(i) +W*w(i)
           W=zeros(2);
           ztt=X2;
           n1=4;
           obj.period= size(X2,2);
           
           % Step1. recover latent variables
           ptt1(:,:,1)= inv(eye(n1)-obj.M)*obj.C*obj.C';
           xtt1(:,1)= zeros(4,1);
           i=1;
           kt(:,:,i)=ptt1(:,:,i)*obj.G'*inv(obj.G*ptt1(:,:,i)*obj.G'+W*W');
           ptt(:,:,i)=ptt1(:,:,i)-ptt1(:,:,i)*obj.G'*inv(obj.G*ptt1(:,:,i)*obj.G'+W*W')*obj.G*ptt1(:,:,i);
           ptt1(:,:,i+1)=obj.M*ptt(:,:,i)*obj.M'+obj.C*obj.C';
           xtt1(:,i+1)= obj.M*xtt1(:,i);
           xtt(:,i)= obj.M*zeros(4,1) +kt(:,:,i)*(ztt(:,i)-obj.G*xtt1(:,i));
           
           for i=2:obj.period
               kt(:,:,i)=ptt1(:,:,i)*obj.G'*inv(obj.G*ptt1(:,:,i)*obj.G'+W*W');
               ptt(:,:,i)=ptt1(:,:,i)-ptt1(:,:,i)*obj.G'*inv(obj.G*ptt1(:,:,i)*obj.G'+W*W')*obj.G*ptt1(:,:,i);
               ptt1(:,:,i+1)=obj.M*ptt(:,:,i)*obj.M'+obj.C*obj.C';
               xtt1(:,i+1)= obj.M*xtt1(:,i);
               xtt(:,i)= obj.M*xtt(:,i-1)+kt(:,:,i)*(ztt(:,i)-obj.G*xtt1(:,i));
           end
           
           output= xtt; % Estimated latent variables
           
           % Step2. Compare the estimated value and realized value
           % Unpack realized values:
               a= X1(1,:);
               z= X1(2,:);
               v= X1(3,:);
               u= X1(4,:);
           % Plot
           time=1:1:obj.period;
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
           hold off
           
           suptitle(name)
           print(str,'-djpeg','-r600')
           snapnow
           
       end
       
       function obj = solve2(obj) % Numerically find the solution
          A0=[1,0,0,0,0,0;
              0,1,0,0,0,0;
              0,0,1,0,0,0;
              0,0,0,1,0,0;
              0,0,0,0,obj.beta,0;
              0,0,0,0,1/obj.sigma,1];
          A1=[obj.rhoa,0,0,0,0,0;
              0,obj.rhoz,0,0,0,0;
              0,0,obj.rhov,0,0,0;
              0,0,0,obj.rhou,0,0;
              0,0,0,1,1,-obj.kappa;
              (1-obj.rhoa/obj.sigma),-(1-obj.rhoz)/obj.sigma,1/obj.sigma,obj.psiy/obj.kappa*obj.sigma,...
              obj.psipi/obj.sigma,1+obj.psiy/obj.sigma];
          C1=[sqrt(diag([obj.sigmaa2, obj.sigmaz2, obj.sigmav2, obj.sigmau2])); zeros(2,4)];
          A=A0\A1;
          obj.C=A0\C1;
          obj.C=obj.C(1:4,1:4);
          
          cutoff =1; % cutoff of stable/unstable variables
          egen= abs(eig(A)) < cutoff;
          n1=4; % number of predetermined variables;
          n2=2; % number of jump variables;
          n= n1+n2; % total number of variables;
          
          % Complex generalized Schur decomposition
          
          [ZZ,TT]=schur(A); % Schur Decompostion to find A=UTU*
          [Z,T] = ordschur(ZZ,TT,'UDI'); % Reorder to put unstable part on the RLS          
          logcon = abs(diag(T)) <= cutoff ;
          
          if sum(logcon) <n1
              warning("Too few stable roots: no stable solutions.");
              obj.M= NaN; obj.G=NaN; J0=NaN;
          elseif sum(logcon) >n1
              warning("Too many stable roots: infinite number of stable solutions.");
              obj.M= NaN; obj.G=NaN; J0=NaN;
          end
          
          Z11=Z(1:n1,1:n1);
          Z21=Z(n1+1:n,1:n1); % Divide Z into a block matrix [Z11; Z21]
          
          obj.M= real(A1(1:n1, 1:n1)); % X1(t+1)= MX1(t)+ e(t+1)
          obj.G= real(Z21/Z11); % X2(t+1)= GX1(t)
       end
   end
end


