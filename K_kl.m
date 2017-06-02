function y = K_kl(a,T)
y = 0.5*(exp((-a*1000)./(8.314*T))); %normalize from kJ/mol to J/mol (?)