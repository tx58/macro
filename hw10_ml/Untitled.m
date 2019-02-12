fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
lb = [0,0.2];
ub = [0.5,0.8];
A = [];
b = [];
Aeq = [];
beq = [];
x0 = [1/4,1/4];
circlecon= @(x)(x(1)-1/3)^2 + (x(2)-1/3)^2 - (1/3)^2;
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,circlecon)



