function clear= mktclear(lamda)

global beta1 beta2 beta3 e1 e2 e3
lamda1=lamda(1);
lamda2=lamda(2);
lamda3=lamda(3);
clear=zeros(3,1);
A=1/(e1+e2+e3)*(lamda1/(1-beta1) + ...
    lamda2/(1-beta2) + lamda3/(1-beta3) );
clear(1)= lamda1/(1-beta1)- e1*A;
clear(2)= lamda2/(1-beta2)- e2*A;
clear(3)= lamda3/(1-beta3)- e3*A;
clear(4)= lamda1+lamda2+lamda3 -1;
end
