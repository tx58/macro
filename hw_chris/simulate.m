function [w, a, x] = simulate(fspace, c, K, w)
global se ro r ygrid yPP
a= zeros(K,1);
%w= zeros(K,1);
x= zeros(K,1);
%w(1)= (length(P)+1)/2;
x(1)= funeval(c,fspace,[0 ygrid(w(1))]);
a(1)= exp(ygrid(w(1)))-x(1);

for t=1:K-1 
    % w(t+1)= ro*w(t)+randn(1)*se;
    % w(t+1)= 1+sum(rand(1)>cumsum(P(w(t),:)));
    x(t+1)= funeval(c,fspace,[(1+r)*a(t) ygrid(w(t+1))]);
    a(t+1)= (1+r)*a(t)+exp(ygrid(w(t+1)))-x(t+1);
end