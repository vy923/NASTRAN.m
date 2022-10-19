function [S,B,P,n] = genFMCparams(A,f)

% COMMENTS BELOW LEFT FROM THE SCRIPT VERSION OF GENFMCPARAMS

% Get bdf file
% A = textin('C:\Users\vvyot\Desktop\NovaSAR\NVBB_01_04_m.dat');

% ------------------ Process input file ------------------
% f = textin('C:\Users\vvyot\Desktop\testpar\test.dat');
% f = textin('C:\Users\vvyot\Desktop\NovaSAR\Model data\FMC_pert_values.dat');

% Remove comments, empty lines and spaces
tmp = contains(f,'%');
f(tmp) = extractBefore(f(tmp),'%');                                         % Removes any '% ...' substrings 
f = strtrim(f);                                                             % Remove leading/trailing whitespace
f = f(~strlength(f)==0);                                                    % Removes lines left empty
    
% Process command continuations
while true
	idx = find(endsWith(f,"..."));
	if any(idx)
        idxm = idx(~endsWith(f(idx+1),"..."));
        f(idxm) = join([extractBefore(f(idxm),"...") f(idxm+1)]);          
        f(idxm+1) = [];
    else
        break
	end
end
    
clearvars idx* tmp

% ---------------------- Read the bdf ------------------------
maxw     = eval(f(1));                                                      % Max written number width for normal/2x width fields. 
cards    = f(3:4:end);
fields   = f(4:4:end);                                                      % As for xbdfread
pertType = f(5:4:end);                                                      % Gaussian, flat
pertVal  = f(6:4:end);                                                      % in terms of CoV/max for normal/flat

n = eval(f(2));                                                             % Number of MC realisations
for k = 1:numel(cards)
    PField{k} = eval(fields{k});
    PType{k}  = eval(pertType{k});
    PVal{k}   = eval(pertVal{k});
    
    [R{k},indA{k},maskQ{k},maskN{k},maskS{k},outf{k}] = xbdfread(A,cards(k),PField{k});
end

idxk = find(~cellfun(@isempty,R));                                          % Index of non-empty cells (e.g. unique cards not missing from the bdf)
npar = numel(idxk);                                                         % Actual number of present unique parameter cards

% --------------------- Clear variables ----------------------
R       = R(idxk);                                                          % Remove empty cells
indA    = indA(idxk);
maskQ   = maskQ(idxk);
maskN   = maskN(idxk);
maskS   = maskS(idxk);
outf    = outf(idxk);
PField  = PField(idxk);
PType   = PType(idxk);
PVal    = PVal(idxk);

clearvars fields pert* f

% --------------- Unperturbed parameter bdfs -----------------
for k = 1:npar                                                              % Counter for cell arrays from xbdfread
    tmp = maskN{k}|maskS{k};
    tmpA{k} = indA{k} + [tmp.*(1:size(tmp,2))-1];                           % Matrix with continuation line indices in A, size = maskS (or maskN)
    tmpA{k} = tmpA{k}.*tmp;                                                 % Format of each row is [l l+1 ... l+h 0 0], where l:l+h are the line indices of the parameter card
    
    tmpAs{k} = nonzeros(tmpA{k}');                                          % Sorted vector of line indices in A containing all k-th parameter data
    C{k} = repmat(padf(A(tmpAs{k}),72),[n 1]);                              % Temporary file w/ n copies of the original k-th parameter subfile, padded to 72 field width
end

B = A;
B(cell2mat(tmpAs')) = [];                                                   % Original bdf with parameter lines removed
P = A(sort(cell2mat(tmpAs')));                                              % Original parameter lines file

% ----------------- Build random parameters ------------------
Rrnd = cell(1,npar);
for k = 1:npar
	Rrnd{k} = distrGen(R{k},PType{k},PVal{k},n);                            % Generate the random parameter values from the required distributions
	RrndC{k} = permute(Rrnd{k},[1 3 2]);
	RrndC{k} = reshape(RrndC{k},size(R{k},1)*n,size(R{k},2));               % Rrnd in format matching C{k}
end

% ---------------- Perturbed parameter bdfs ------------------
% Set field indices
f8  = {1:8 9:16 17:24 25:32 33:40 41:48 49:56 57:64 65:72};
f16 = {1:8 9:24 25:40 41:56 57:72 73:88 89:104 105:120 121:136};  

% Loop over cards
for k = 1:npar

	maskQs  = sum(maskQ{k},2);                                              % Total number of lines each parameter card takes
	maskQcs = cumsum(maskQ{k},2);                                           % Offset from first to l-th line of the cards (due to 8/16 field format)
    maskQcsr = repmat(maskQcs,[n 1]);                                       % The above mask copied n times
    
	summQ2 = sum(maskQ{k},2);                                               % Total number of lines taken by each card 
	idxl = cumsum(summQ2);                                                  % Index of the last line of the i-th card in A(tmpAs{k})
	idxl = idxl - summQ2 + 1;                                               % Index of the first line -//-
    
	totl = sum(summQ2);                                                     % Total number of lines of A(tmpAs{k})
	idxlr = idxl + [0:totl:totl*(n-1)];                                     % First line index in C{k}
	idxlr = idxlr(:);                                                       % -//- expanded as a column vector
    
    % Loop over card lines
	cc = 0;                                                                 % Counter for columns in Rrnd
	for l = 1:min(numel(outf{k}),size(maskQ{k},2))                          % Line of the param card / UPDATE [14/5/18]: was numel(outf{k}) only
        
        mask1f = maskQ{k}(:,l)==1;                                          % Fields in the l-th line in 8-char format
        mask2f = maskQ{k}(:,l)==2;                                          % -//- 16-char format
        
        chk(1) = any(mask1f);
        chk(2) = any(mask2f);
 
        mask1r = repmat(mask1f,[n 1]);                                      % 8 character lines in RrndC{k}
        mask2r = repmat(mask2f,[n 1]);                                      % 16 character lines in RrndC{k}

        % Loop over fields in each line
        for t = 1:numel(outf{k}{l})                                         % Loop over requested fields of the l-th line
            cc = cc + 1;                                                    % R-column counter increment
            
            if chk(1)
                fval1  = RrndC{k}(mask1r,cc);                               % cc-th column of 8-fld values
                fval1s = padf(num2field(fval1,maxw(1)),8);                  % Double-to-char formatted conversion - see num2field/padf core functions
                if l>=2     
                    idxla = idxlr(mask1r) + maskQcsr(mask1r,l-1);           % Index adjustment for updating card fields on line number >= 2
                else
                    idxla = idxlr(mask1r);
                end
                C{k}(idxla,f8{outf{k}{l}(t)}) = fval1s;                     % Update C{k} with the new parameter values
            end 
            
            if chk(2)                                                       % Same procedure as for normally sized fields
                fval2  = RrndC{k}(mask2r,cc);
                fval2s = padf(num2field(fval2,maxw(2)),16);
                if l>=2     
                    idxla = idxlr(mask2r) + maskQcsr(mask2r,l-1);
                else
                    idxla = idxlr(mask2r);
                end   
                if outf{k}{l}(t) > 5                                        % Additional index adjustment for fields in the second part of a *-line
                    C{k}(idxla+1,f16{outf{k}{l}(t)-4}) = fval2s;
                else
                    C{k}(idxla,f16{outf{k}{l}(t)}) = fval2s;
                end
            end
            
        end % FOR fields
    end % FOR lines

   % Reshape each block in C back to individual files (small overhead)
   S{k} = reshape(C{k}',size(C{k},2),totl,n);
   S{k} = permute(S{k},[2 1 3]);
   
end % FOR parameters

% Assemble S{k} into final parameter bdfs ready for output
S = vertcat(S{:}); 

end % function genFMCparams

% ------------------------ Functions -------------------------
% Computes an array of random parameter values 
function Rrnd = distrGen(R,type,cov,N)
    if nargin<=3 
        N = 1;
    end

	maskG = type==1;
    maskF = type==2;

    nrow = size(R,1);                                                       % Rows (i.e. number of times the card is repeated in the bdf)
    nG = sum(maskG);                                                        % Number of fields in the card with normal distribution (i.e. E,nu,...)
    nF = sum(maskF);                                                        % Number of fields with flat distribution 
    
    distG = normrnd(ones(nrow,nG,N), ones(nrow,nG,N).*cov(maskG));
    Rrnd(:,maskG,:) = distG.*R(:,maskG);                                    % Normally distributed values 
    
    distF = 1 + (1-2*rand(nrow,nF,N)).*cov(maskF);
    Rrnd(:,maskF,:) = distF.*R(:,maskF);                                    % Flat distribution values 
end




%{
% Output all files with the textout routine. 'deblank'=true is optional,
% but recommended, as it reduces output file size
for i = 1:size(S,3)
    textout(S(:,:,i),['C:\Users\vvyot\Desktop\testpar\params\p' num2str(i) '.dat'],true);
end

% Simple method to find all unique cards in A
    Ac = A(~startsWith(A,{char(0),' ','*','$','+'}));
    Ac = Ac(Ac~="");
    Ac = char(Ac);
    Ac = Ac(:,1:8);
    Ac = unique(Ac,'rows');

% ------- EXAMPLE INPUT FILE -------
% Text formatting below is intentionally poor, for testing the above code

% Max number characters to for 8/16 char fields; [7,15] makes the ourput file 
% clearer, but [8,15/16] is recommended for best number representation accuracy.
	[7  15]

% Number of instances
	2e3
	
%
% Set perturbation parameters: nearly native ML format is acceptable
MAT1
[3:6]               	% E, G, nu, rho
ones(1,4)				% Gaussian distribution
[0.07 0.05 ...			% Coefficient of variation

...

 0.02 ...
 0.0] 	
 
%
PBUSH
[4:9]               	% K1:K6
2*ones(1,6)				% Flat distribution
[zeros(1,3) .2 .2 .2]	% No perturbation for first 3 values
%
PSHELL
[3,9]               	% T, NSM
[1 2]					% [Gaussian Flat]
[.1 .1]
%
% MAT9 not contained in the test bdf, intentionally put here for 
% verification if the code works when missing cards are requested
%
MAT9					
{[3:9];[2:9];[2:8]}	 	% Gij(i=1:6,j=1:6,i>=j), rho
ones(1,22)				% Gaussian
[0.05*ones(1,21) 0.02]
%
CONM2
{[5];[2:7]}         	% M; Iij(ij = 11,21,22,31,32,33)
ones(1,7)				% Gaussian
[0.03 0.08*ones(1,6)]
%
%
CBAR					
{[6:8];[2:9]}	 		% Some values w/ no reasonable meaning, only added here for testing
ones(1,11)				% Gaussian
[0.05*ones(1,9) .02 .1] % Arbitrarily selected pert. values, for testing

%}
