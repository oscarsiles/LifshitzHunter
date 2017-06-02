function n3 = calcN(n1,n2,f1)
%CALCN Calculate n for all binary mixures

% boundary conditions
if f1 == 1.0
	n3 = n1;
	return
elseif f1 == 0.0
	n3 = n2;
	return
end

if n1 == n2 % if the same liquid, refractive index won't change
	n3 = n1;
	return
else
	f2 = 1 - f1;
	% Lorentz-Lorenz
	%n3 = abs((sqrt(-2.*f1.*(n1^2-1).*(n2^2+2)-(n1^2+2).*(2.*f2.*(n2^2-1)+n2^2+2)))./(sqrt(f1.*(n1^2-1).*(n2^2+2)+(n1^2+2).*(f2.*(n2^2-1)-n2^2-2))));
	% Heller's relation (faster, equivalent results)
	m = n2./n1;
	n3 = n1.*(1.5.*f2.*((m^2-1)./(m^2+2)))+n1;
end

end

