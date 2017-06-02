function [Wad3,Wad5] = calcWAD(inputSurfaces,surface1,surface2,eSolvent,nSolvent,T)
%CALCWAD Calculate work of adhesion for Lifshitz theory.
D_0 = 0.165*10^(-9); % from israelachvilli

% import parameters from surface
surface1Row = ismember(inputSurfaces,surface1);
surface2Row = ismember(inputSurfaces,surface2);

e1 = sscanf(cell2mat(inputSurfaces(surface1Row,2)), '%f');
e2 = sscanf(cell2mat(inputSurfaces(surface2Row,2)), '%f');
n1 = sscanf(cell2mat(inputSurfaces(surface1Row,3)), '%f');
n2 = sscanf(cell2mat(inputSurfaces(surface2Row,3)), '%f');

% check for "empty"
if isnan(e1) || isnan(e2) || isnan(eSolvent) || isnan(n1) || isnan(n2) || isnan(nSolvent)
	Wad3 = '-';
    Wad5 = '-';
	return
end

% 3-medium Hamaker and calculation
A132 = calcA(e1,e2,eSolvent,n1,n2,nSolvent,T);
Wad3 = A132./(12*pi.*D_0^2);


% import parameters for 5-medium calculation
eUndecane = sscanf(cell2mat(inputSurfaces(surface1Row,7)), '%f');
nUndecane = sscanf(cell2mat(inputSurfaces(surface1Row,8)), '%f');

molVol1 = sscanf(cell2mat(inputSurfaces(surface1Row,9)), '%f');
eEndGroup1 = sscanf(cell2mat(inputSurfaces(surface1Row,10)), '%f');
nEndGroup1 = sscanf(cell2mat(inputSurfaces(surface1Row,11)), '%f');

molVol2 = sscanf(cell2mat(inputSurfaces(surface2Row,9)), '%f');
eEndGroup2 = sscanf(cell2mat(inputSurfaces(surface2Row,10)), '%f');
nEndGroup2 = sscanf(cell2mat(inputSurfaces(surface2Row,11)), '%f');

if isnan(molVol1) || isnan(eEndGroup1) || isnan(nEndGroup1) || isnan(molVol2) || isnan(eEndGroup2) || isnan(nEndGroup2) || isnan(eUndecane) || isnan(nUndecane)
    Wad5 = '-';
    return
end

% calculate thickness of end group layer (assume sphere, double the radius)
T1 = nthroot((molVol1.*(3/4)./(pi.*6.022E23)),3).*2;
T2 = nthroot((molVol2.*(3/4)./(pi.*6.022E23)),3).*2;

% find all Hamaker constants for force calc (5-medium)
A232d = calcA(eEndGroup1,eEndGroup2,eSolvent,nEndGroup1,nEndGroup2,nSolvent,T);
A121 = calcA(eUndecane,eUndecane,eEndGroup1,nUndecane,nUndecane,nEndGroup1,T);
A32d3 = calcA(eSolvent,eSolvent,eEndGroup2,nSolvent,nSolvent,nEndGroup2,T);
A1d2d1d = calcA(eUndecane,eUndecane,eEndGroup2,nUndecane,nUndecane,nEndGroup2,T);
A323 = calcA(eSolvent,eSolvent,eEndGroup1,nSolvent,nSolvent,nEndGroup1,T);

% 5-medium calculation
fun = @(x) ((1/(6.*pi)).*(A232d./(x.^3)-(sqrt(A121.*A32d3))./((x+T1).^3)-(sqrt(A1d2d1d.*A323))./((x+T2).^3)+(sqrt(A1d2d1d.*A121))./((x+T1+T2).^3))); % x = D
Wad5 = integral(fun,D_0,Inf); % check limits?

end