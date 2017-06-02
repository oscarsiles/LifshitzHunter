function y = E_kl(a,b,T)
y = (a.*b)./(1+exp((a.*b*1000)./(8.314*T)))-5.6; % kJ/mol