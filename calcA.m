function A = calcA(e1,e2,e3,n1,n2,n3,T)
%CALCA Calculate Hamaker constant for 3-medium system.
v_e=3.3*10^15;

if e1 == e2 && n1 == n2 % if same surfaces, use shorter equation (cheaper)
	A = ((3*1.38*10^(-23).*T)/4).*((e1-e3)./(e1+e3))^2+((3*6.63*10^(-34).*v_e)/(16*sqrt(2))).*((n1^2-n3^2)^2/(n1^2+n3^2)^(3/2));
else
	A = ((3*1.38*10^(-23).*T)/4).*((e1-e3)./(e1+e3)).*((e2-e3)./(e2+e3))+((3*6.63*10^(-34).*v_e)/(8*sqrt(2))).*(((n1^2-n3^2).*(n2^2-n3^2))./((sqrt(n1^2+n3^2).*sqrt(n2^2+n3^2)).*(sqrt(n1^2+n3^2)+sqrt(n2^2+n3^2))));
end

end