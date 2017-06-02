function dG_kJ = calcDG(filename, K, T)
%CALCDG Calculates dG for Hunter model. Output in kJ/mol
if exist(filename, 'file') == 0 %if file does not exist
	dG_kJ = '-';
	return
end

c_max = 300;
phaset = read_mixed_csv(filename, ',');

% values of dG_s will be in string format
dG_si = str2double(phaset(224,19))*1000; %convert string to number, x1000 to convert into J from kJ
dG_sj = str2double(phaset(225,19))*1000;
theta = str2double(phaset(224,17));

% calculate dG_c and dG_b
dG_c = -2*8.314.*T.*log((sqrt(1+8*theta)-1)./(4*theta));
dG_b = 2*8.314.*T.*log((sqrt(1+8*K.*theta)-1)./(4*K.*theta));

% put it all together
dG = dG_b + dG_c - dG_si - dG_sj + 8.314.*T.*log(c_max);
dG_kJ = dG/1000;
end

