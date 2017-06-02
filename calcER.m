function e3 = calcER(e1,e2,f1)
%CALCER Calculate e_r for binary mixtures

% boundary conditions
if f1 == 1.0
	e3 = e1;
	return
elseif f1 == 0.0
	e3 = e2;
	return
end

if e1 == e2
	e3 = e1;
	return
end

e3 = f1.*e1+(1-f1).*e2;

% if p1 < 0.2 && p2 < 0.2 % non-polar mixture, use linear mixing
% 	e3 = f1.*e1+(1-f1).*e2;
% 	return
% else % either polar-nonpolar or polar-polar
% 	e3 = f1.*e1+(1-f1).*e2; % placeholder for now
% 	return
% end

end

