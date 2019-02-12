function T=DSGE_solve(tau,beta,theta,phi_pi,phi_y,varphi,alpha,eps,rho_v,rho_a,rho_z,rho_u,sigma_v,sigma_a,sigma_z,sigma_u)

    % Unpack parameters
    sigma= tau; % consumption elasiticity
    phi= varphi; % labor elasticity
    epsilon= eps; % different goods bundle elasticity
    alpha= alpha; % production coefficient
    beta= beta; % intertemporal discount factor
    theta= theta; % ratio of firms which cannot change prices
    psipi= phi_pi; % taylor rule of inflation
    psiy= phi_y; % taylor rule of output
    rhoa= rho_a; % technology AR(1) coefficient
    sigmaa2= sigma_a^2; % technology AR(1) variance
    rhoz= rho_z; % demand shock AR(1) coeffeicient
    sigmaz2= sigma_z^2; % demand shock AR(1) variance
    rhov= rho_v; % policy shock AR(1) coefficient
    sigmav2= sigma_v^2; % policy shock AR(1) variance
    rhou= rho_u; % cost-push shock AR(1) coefficient
    sigmau2= sigma_u^2; % cost-push shock AR(1) coefficient

    rho= -log(beta);  % linearied beta
    psiya= (1+phi)/(sigma*(1-alpha)+phi+alpha); % tech effect on output, affecting efficient level of output 
    lamda= (1-theta)*(1-beta*theta)/theta*(1-alpha)/(1-alpha+alpha*epsilon);
    kappa= lamda*(sigma+(phi+alpha)/(1-alpha));
    gamma= 1/((1-beta*rhoa)*(sigma*(1-rhoa)+psiy)+ kappa*(psipi-rhoa));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Stable/unstable Uncoulping to solve the model
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
    T= real(Z21/Z11); % X2(t+1)= GX1(t) 


end