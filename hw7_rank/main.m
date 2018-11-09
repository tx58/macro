modela=NKM;
modela=init(modela);

% introduce a one-time shock:
v=zeros(20,1);
modela.shock(v, 'q4_1', 'Steady state without positive technology shocks');
v(2)=0.1;
modela.shock(v, 'q4_2', 'Impulse Response for a positive technology shock: Regime 1');
modelb=NKM;
modelb.psipi=10;
modelb.psiy=0;
modelb=init(modelb);
modelb.shock(v, 'q6_1', 'Impulse Response for a positive technology shock: Regime 2');

% introduce a 100 period random shock:
v= sqrt(modela.sigmaa2)*randn(100); % tech shock
modela.shock(v, 'q5_1', 'Path for 100 period random shocks: Regime 1');
modelb.shock(v, 'q5_2', 'Path for 100 period random shocks: Regime 2');