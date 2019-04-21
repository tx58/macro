%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negeshi method when \beta is small 
% Author: Tianli Xia
% Feb 17th, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem setup:
global beta1 beta2 beta3 e1 e2 e3
beta2=0.95
beta3=0.95
e1=0.1
e2=1
e3=2
t=20

for beta1= 0.05
    lamda = fsolve( @mktclear, [1 1 1] )
    beta= [beta1 beta2 beta3];
    c = lamda'.*beta'.^repmat(1:1:t,3,1)./ ...
        sum(lamda'.*beta'.^repmat(1:1:t,3,1)).*(e1+e2+e3);
end

plot(1:1:t, c(1,:), "-k")
hold on
plot(1:1:t, c(2,:), '--k')
hold on
plot(1:1:t, c(3,:), '-.k')
hold on
legend("\beta1=0.05,e1=0.1", "\beta2=0.95,e2=1", "\beta3=0.95,e3=2", 'location', 'best')
xlabel("period")
ylabel("consumption")
title("Income distribution")
print -djpeg -r600 hw4_income_distribution
