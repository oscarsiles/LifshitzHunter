function combined(surface1,surface2)
% Calculates dG (Hunter) and Wad (Lifshitz) for two surfaces in solvent(s)
if nargin == 0 % no input arguments
    disp('No input arguments, defaulting to dodecane/dodecane.');
    surface1 = 'dodecane';
    surface2 = 'dodecane';
elseif nargin == 1
    disp('Wrong number of input arguments.');
    return
end

% Declare constants
T=298;

% import data (use custom function to create a mixed cell array)
input = read_mixed_csv('./input/inputParameters.csv', ',');
inputSurfaces = read_mixed_csv('./input/surfaces.csv', ',');
surface1Row = ismember(inputSurfaces,surface1);
surface2Row = ismember(inputSurfaces,surface2);

outputFolder = strcat(pwd,'/output/',surface1,'/',surface2,'/');

% check that output folders exist, if not create it
if exist(outputFolder, 'dir') ~= 7 % 7 is type=folder
	mkdir(outputFolder)
end

% copy data before manipulation and find some information about input size
outputSingle = input;
[rows,columns] = size(input);

% Lifshitz import data
e_i = sscanf(cell2mat(inputSurfaces(surface1Row,4)), '%f'); % look into input arguments for these (review has values for donor (+, 2.7 for alcohol) and acceptor (-))
e_j = sscanf(cell2mat(inputSurfaces(surface2Row,5)), '%f');

% Hunter preliminary calculations
% find all combinations of E/K for i and j (i and j interchangeable)
for i = 1:3
    if i == 1 %E_ij
        E_ij = E_kl(e_i,e_j,T);
        K_ij = K_kl(E_ij,T);
    elseif i == 2 %E_ii
        E_ii = E2_kl(e_i,e_i);
        K_ii = K_kl(E_ii,T);
    elseif i == 3 %E_jj
        E_jj = E2_kl(e_j,e_j);
        K_jj = K_kl(E_jj,T);
    end
end
% find total K (K_ij = K_ji)
K = (2*K_ij+K_ii+K_jj)/4;


% Single solvent

% write column labels (single)
outputSingle{1,columns+1} = 'Wad';
outputSingle{1,columns+2} = 'Wad (5-medium)';
outputSingle{1,columns+3} = 'dG';

% loop for all solvents (start at row 2, row 1 is labels)
for j = 2:rows
    %Lifshitz
    eSolvent = str2double(input(j,3));
    nSolvent = str2double(input(j,4));
	[Wad3,Wad5] = calcWAD(inputSurfaces,surface1,surface2,eSolvent,nSolvent,T); % calculate work of adhesion for 3- and 5-medium systems
    outputSingle{j,columns+1} = Wad3;
    outputSingle{j,columns+2} = Wad5;
    
    % Hunter
    % read solvent name and import corresponding csv from phaset subfolder
    solvent1 = char(input(j,2));
	fprintf(strcat(solvent1,'..'));

	filename = strcat(pwd,'/output_pt/',surface1,'/',surface2,'/single/phaset_',solvent1,'.csv');
	
	dG_kJ = calcDG(filename, K, T); % output in kJ
    outputSingle{j,columns+3} = dG_kJ;
    fprintf('Done.\n');
end
% export single solvent
cell2csv(strcat(outputFolder,'outputSingle.csv'), outputSingle, ',');


% Mixtures

% mixture constants
n = 9; % number of fractions (9 for 0.1-0.9 in 0.1 increments)
i = 2; % start counter (1 is label)
counter = 1; % start counter, to see if all mixtures calculated
% write column labels (mixture)
outputMix = cell(1);
outputMix{1,1} = strcat('e_i=',num2str(e_i,2)); % 2 significant figures
outputMix{1,2} = strcat('e_j=',num2str(e_j,2));
outputMix{1,4} = 'Wad';
outputMix{1,4+n+i} = 'dG';
outputMix{2,1} = 'solvent 1+2';
% outputMix{2,1} = 'solvent 1';
% outputMix{2,2} = 'solvent 2';

% loop for all possible binary mixtures, from 0.0 to 1.0 of FIRST solvent listed in 0.1 steps
for k = 2:(rows-1) % start from first solvent, until second to last solvent (line 1 is header)
    for l = 1:(rows-k) % combine all solvents, stopping at last one
		% read parameters
        solvent1 = char(input(k,2));
        solvent2 = char(input(k+l,2));
		eSolvent1 = str2double(input(k,3));
		eSolvent2 = str2double(input(k+l,3));
		nSolvent1 = str2double(input(k,4));
		nSolvent2 = str2double(input(k+l,4));
% 		pSolvent1 = str2double(input(k,4));
% 		pSolvent2 = str2double(input(k+l,4));
		
		disp(strcat(solvent1,'+',solvent2));
        outputMix{i+1,1} = strcat(solvent1,'+',solvent2);
		%outputMix{i+1,1} = solvent1;
        %outputMix{i+1,2} = solvent2;
        for m = 0:(n+1) % include the volume fractions				
			% create correct fraction for input
			if m == (n+1)
				fraction = '1.0';
				f1 = 1.0;
			else
				fraction = strcat('0.',int2str(m));
				f1 = m/(n+1.0);
			end
			fprintf(strcat(fraction,'..'));
			
			% Lifshitz
% 			eMix = calcER(eSolvent1,eSolvent2,pSolvent1,pSolvent2,f1);
            eMix = calcER(eSolvent1,eSolvent2,f1);
			nMix = calcN(nSolvent1,nSolvent2,f1);
			[Wad3Mix,Wad5Mix] = calcWAD(inputSurfaces,surface1,surface2,eMix,nMix,T);
			%outputMix{i+1,columns+m} = Wad3Mix;
            outputMix{i+1,columns+m} = Wad5Mix; % use 5-medium. error within 5%!
			
			% write column labels (+1 because starting at m=0) on first pass
			if (i == 2)
				outputMix{2,columns+m} = fraction;
				outputMix{2,columns+m+n+2} = fraction;
			end
			
			filename = strcat(pwd,'/output_pt/',surface1,'/',surface2,'/mixtures/',solvent1,'/',solvent2,'/phaset_',solvent1,'_',solvent2,'_',fraction,'.csv');
			dG_kJ = calcDG(filename, K, T);
            outputMix{i+1,columns+m+n+2} = dG_kJ;
        end
		fprintf('Done.\n');
        i = i+1; 
		counter = counter + 1; % for debugging (check if right number of input files)
    end
end

% export binary mixture
cell2csv(strcat(outputFolder,'outputMix.csv'), outputMix, ',');
disp(strcat(int2str(counter),' combinations processed.'));
