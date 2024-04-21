function [R,indA,maskQ,maskN,maskS,outf] = xbdfread(A,card,outFields)

%{
    > Finds lines containing 'card' in a text file
    > Finds corresponding row indices and continuation lines
    > The card data requested by outFields is read onto a matrix, or cell
	array in the case of RBE2
    > Fast Nastran num format (no exp. sign) -> double by str2val function
    > Returns R in the same order as the cards are found in A
    > Adding new cards only requires an extra 'case' in the cardfields function
    
    > outFields can be 
        (1) a cell array, e.g. outFields{i} = 2:5 would return the 2,3,4,5th
            entries of the i-th line of each matched card
        (2) a logical matrix of size max. (lines<=max card lines) x (fields<=8)
            e.g. for GRID, max card lines = 1, 3 for CHEXA, etc.
    > For entries of nonconsecutive lines: outFields = {[2:4];[];[5 7]}
    > indA has the indices of rows containing 'card'
    > See code below for maskQ, maskN, maskS

    - RBE2 can only be read from 8-character fields in the current version
    - PBUSH: Assumes first line is always stiffness (flag 'K')
    - PBEAML: Only reads first line (unlimited length card, like the RBEs)
    - PBARL: Same comment as for PBEAML
%}

% Prevent getting error from bdfscan if last file line is nonempty
	if isstring(A),     A(end+1) = "";
	elseif iscell(A),   A{end+1} = '';
    end
    
% Convert to char is string or cell
    if ~ischar(card)    
        card = char(card);    
    end
    
% Covert outFields to the necessary cell format
    try
        if ~iscell(outFields)                                               % Convert from matrix format to cell of indices, if it is a logical array
            isbool = all(ismember(outFields(:),[0,1]));
            for i = 1:size(outFields,1)
                if isbool,  tmp{i} = find(outFields(i,:));                  % Each row of logical values corresponds to a line of the card
                else,       tmp{i} = outFields(i,:);
                end
            end
            outFields = tmp;
        end   
    catch
        outFields = [];                                                     % If the var does not exist, default outputs are used 
    end

% Set default values for outFields if the read text field is empty
	[empt,outf] = cardfields(card,outFields);

% Search file - find line indices of the required card
    maskA = startsWith(A,card);
    indA  = find(maskA);

% Extract continuation line indices, if any matching cards are found 
    if isempty(indA)
        emptyVars = cell(6,1);
        [R,indA,maskQ,maskN,maskS,outf] = emptyVars{:};
     %   R = [];
     %   maskQ = [];
        return                                                              % Break execution
    else
        [maskNS,maskS] = xbdfscan(A,indA,{' ','+','*'},card);               % maskS = true for long field lines, maskNS for either
     	maskN = maskNS & ~maskS;                                            % maskN = true only for short field lines of 'card'
    end

% ----- Separate computations for RBE2 vs other cards -----   
    if ~any(contains({'RBE2','RBE3'},card))

    % Initialise R to an array of size indA by number of requested fields
        emf = cell2mat(empt');                                              % All empty field values in a vector
        n   = length(indA);                                                 % Total number of matched cards
        R   = zeros(n,numel(emf)) + emf;                                    % Set to default card values

    % Compute a mask of 0,1,2s, assuming (*) lines are always in pairs
        maskQ   = zeros(size(maskN));                                       % Initially let maskQ have the maximum possible size
        posQ    = (1:n)';                                                   % Linear index vector of the next positions in maskQ to update
        tmpS    = maskS(:,1);                                               % Ensures [1 1]->[2 0], but [1 1 1]\->[2 2 0] from maskS to maskQ

        for i = 1:size(maskN,2)  
            mNi = maskN(:,i);                                               % For faster memory access
            maskQ(posQ(mNi)) = 1;                                           % Normal cont. lines always directly go to the next corresponding posQ
            posQ(mNi) = posQ(mNi)+n;                                        % Update posQ

            if i>=2                                                         % (*)-lines
                tmpSS = maskS(:,i-1) & maskS(:,i) & tmpS;                   % Checks for two consecutive 1s in maskS, and no 2 at the last corresponding position in maskQ
                maskQ(posQ(tmpSS)) = 2;
                posQ(tmpSS) = posQ(tmpSS)+n;          
                tmpS(tmpSS) = false;           
                tmpS(~tmpSS & maskS(:,i)) = true;          
            end
        end
        
        maskQ(:,~any(maskQ,1)) = [];                                        % Delete empty columns of maskQ

    % Break if long field lines are not always paired    
        if sum(maskNS(:))~=sum(maskQ(:))
            error('|>| Non-paired (*) lines present |<|')
        end
        
    % Set field indices
        f8  = {1:8 9:16 17:24 25:32 33:40 41:48 49:56 57:64 65:72};
        f16 = {1:8 9:24 25:40 41:56 57:72 73:88 89:104 105:120 121:136};        

    % Read the card data 
        colR = 0;                                                           % Initialise column count for the output R
        csumQ = cumsum(maskQ,2);                                            % Accummulated sum of maskQ, used to compute index offset for indA
            
        for line = 1:numel(outf)                                            % Loop over the requested card lines
        
            lfs = outf{line};                                               % field indices for the current iteration                                                  
            if isempty(lfs) || size(maskQ,2)<line                           % Skip to next iteration if no field of 'line' is requried
                continue
            end
        
            tmpN = maskQ(:,line)==1;                                        % mask used for row indexing of R (non-* lines)
            tmpS = maskQ(:,line)==2;
            atmpN = any(tmpN);                                              % Save computation
            atmpS = any(tmpS);
            
            if atmpN
                offN = csumQ(tmpN,line)-1;                                  % Offset from indA
                AtmpN = padf(A(indA(tmpN)+offN),72);                        % Read the corresponding A lines and expand to 72 width  
            end            
            if any(tmpS)
                offS = csumQ(tmpS,line)-2;
                AtmpS1 = padf(A(indA(tmpS)+offS),72);                       % 1st part of the long field line
                AtmpS2 = padf(A(indA(tmpS)+offS+1),72);                     % 2nd part 
            end
            
            for pos = 1:numel(lfs)
                colR = colR+1;                                                      % Column index of R   
                if atmpN                                                            % Read from normal continuation lines
                    R(tmpN,colR) = str2val(AtmpN(:,f8{lfs(pos)}),emf(colR));        
                end
                if atmpS && lfs(pos)<=5                                             % 1st part of long lines
                    R(tmpS,colR) = str2val(AtmpS1(:,f16{lfs(pos)}),emf(colR));      
                elseif atmpS && lfs(pos)>5                                          % 2nd part
                    R(tmpS,colR) = str2val(AtmpS2(:,f16{lfs(pos)-4}),emf(colR));	% Updated from lfs(pos-4) [3/3/2018]  
                end                 
            end
           
        end % for line=1:numel(outf)
        
% ----- Prepare RBE2 for reading. Only works with 8-char wide fields -----
    else 
        
    % Allocate to avoid output errors
        maskQ = missing;                                                    % MOD 15.02.23 / bdfread -> xbdfread in NGBA03 reductions script

    % Verify no (*)-lines are present
        if any(maskS(:))
            error('|>| Long field lines present in RBE2, syntax not implemented |<|')    
        end
   
	% Group RBE cards by number of continuation lines
        m  = sum(maskN,2);
        mu = unique(m);
        mm = bsxfun(@plus,indA(maskN(:,1)),0:size(maskN,2)-1).*maskN; 
    
    % Loop through the RBE sets by number of continuation lines
        for q = 1:numel(mu)
            indm = mm(m==mu(q),1:mu(q));                                    
            datm = A(indm);
            indN{q} = indm; 
            
            datm    = padf(datm,72);                                        % Convert to char and expand/truncate to 72 fields
            datm    = datm(:,9:size(datm,2),:);
            datm    = reshape(datm,size(datm,1),size(datm,2)*size(datm,3));
            datN{q} = reshape(datm',8,numel(datm)/8)';
        end  

	% Read the RBE2 data onto a cell array
        rbeNum = 1;                                                         % Initialise RBE2 card count
        for k = 1:numel(datN)                                               % Loop through the RBEs as read by number of continuation lines
            tmp = datN{k};
            maskNonGrid = contains(string(tmp),{'.','-','+'});              % A bit slow due to string conversion

            if any(maskNonGrid)
                tmp(maskNonGrid,:) = ' ';                                   % Removes coefficient of thermal expansion data
            end    

            tmp = str2val(tmp);                                             % Converts empty spaces to NaN
            tmp = reshape(tmp,(mu(k)*8),(numel(tmp)/mu(k)/8));              % Each column corresponds to a separate RBE2

            for s = 1:size(tmp,2)                                           % Parse and collect RBE2 data onto a cell array
                tmpG = tmp(4:end,s);
                R(rbeNum,:) = {tmp(1,s) tmp(2,s) tmpG(~isnan(tmpG)) tmp(3,s)};   % {[ID] [master] [dependent] [components]}
                rbeNum = rbeNum+1;
            end  
        end % Loop through datN entries  

    end % if ~any(contains({'RBE2','RBE3'},card))

end % xbdfread 



% ----------------------------- FUNCTIONS ---------------------------

% Default values to read for all accepted cards
function [emptyVals,outFields] = cardfields(card,outFields)
     
% Output data - comments refer to the 'default' fields only
switch card
    case 'CBAR'                                                             % ID, PID, G1, G2
        default = {2:5};
        emptyField = NaN(2,9);    
	case 'CBEAM'                                                            % ID, PID, G1, G2
        default = {2:5};
        emptyField = NaN(3,9); 
    case 'CBUSH'                                                            % ID, PID, G1, G2, CID
        default = {[2:5 9]};
        emptyField = NaN(2,9);
    case 'CHEXA'                                                            % ID, PID, G1:G8
        default = {2:9; 2:3};
        emptyField = NaN(3,9);
    case 'CONM2'                                                            % ID, G1, CID, M, Iij(ij = 11,21,22,31,32,33)
        default = {2:5; 2:7};
        emptyField = NaN(2,9);
    case 'CORD2R'                                                           % CID, CID2, A1:A3, B1:B3, C1:C3
        default = {2:9; 2:4};
        emptyField = [NaN NaN 1 NaN(1,6); NaN(1,9)];        
    case 'CQUAD4'                                                           % ID, PID, G1, G2, G3, G4
        default = {2:7};
        emptyField = NaN(2,9);     
    case 'CROD'                                                             % ID, PID, G1, G2
        default = {2:5};
        emptyField = NaN(1,9);         
    case 'CTRIA3'                                                           % ID, PID, G1, G2, G3
        default = {2:6};
        emptyField = NaN(2,9);    
    case 'GRID'                                                             % ID, CID, X, Y, Z
        default	= {[2 3:6]};
        emptyField = [NaN NaN 1 NaN(1,6)];
    case 'MAT1'                                                             % ID, E, G, nu, rho
        default = {2:6};                                                  
        emptyField = NaN(2,9);
    case 'MAT8'                                                             % ID, E1, E2, nu12, G12, G1z, G2z, rho
        default = {2:9};                                                  
        emptyField = NaN(3,9);
    case 'MAT9'                                                             % ID, Gij(i=1:6,j=1:6,i>=j), rho
        default = {2:9; 2:9; 2:8};                                      
        emptyField = NaN(4,9);        
    case 'PBARL'                                                            % ID, matID, type // dim1, dim2, ..., NSM
        default = {2:3; 2:3};
        emptyField = NaN(2,9);    
    case 'PBEAML'                                                           % ID, matID // dim1, dim2
        default = {2:3; 2:3};
        emptyField = NaN(2,9);         
    case 'PBUSH'                                                            % ID, K1:K6
        default = {[2 4:9]};                                                
        emptyField = NaN(4,9);
    case 'PROD'                                                             % ID, matID, Area, J, C, NSM
        default = {2:7};
        emptyField = NaN(1,9);          
    case 'PSHELL'                                                           % ID, matID1, t, matID2, 12I/T3, matID3, TS/T, NSM
        default = {2:9};
        emptyField = NaN(2,9);
    case 'PSOLID'                                                           % ID, matID, cordM, IN, STRESS, ISOP
        default = {2:7};
        emptyField = NaN(1,9);
    case 'RBE2'
        outFields = [];
        emptyVals = [];
        return 
    otherwise
        error(['|>| Card ' card ' not implemented in xbdfread |<|'])
end % SWICTH card

% Set output to default if custom is not required
    try     outFields(1);
    catch,  outFields = default;
    end   
    
    % Get the default empty values corresponding to outFields
    if length(outFields)>size(emptyField,1)
        error('|>| More output lines requested than the max number of card lines |<|')
    else
        for i = 1:numel(outFields)
            emptyVals{i,1} = emptyField(i,outFields{i});
        end    
    end
end

% Finds indices of continuation lines  
function [maskC,maskS] = xbdfscan(A,indA,flag,card)
    maskC = true&indA;                                                      % Same as true(size(indA))
    maskS = startsWith(A(indA),[card '*']);                                 % Initialise the long field line mask
    
    k = 1;
    while true       
        ind = indA(maskC(:,k))+k; 
        Aind = A(ind);                                                      % To avoid computing A(ind) twice
        varC = startsWith(Aind,flag);                                       % Both 16- and 8-field continuation lines
        varS = startsWith(Aind,'*');                                        % Only 16-field continuation lines
        if any(varC)
            maskC(maskC(:,k),k+1) = varC;
            maskS(maskC(:,k),k+1) = varS;          
            k = k+1;
        else
         	break
        end
    end    
end   

% Faster than pad, can also truncate the char array
function str = padf(str,w)
    if ~exist('w','var'),       w = 72;             end                     % Set desired char array width
    if ~ischar(str),            str = char(str);    end                     % Convert to char array

	if size(str,2)==w                                                       % Do nothing
        return
    elseif size(str,2)<w                                                    % Expand to 72 fields per line, if less
        if      ismatrix(str), 	str(:,size(str,2)+1:w) = ' ';               % Faster than (:,~,:)
        elseif  ndims(str)==3, 	str(:,size(str,2)+1:w,:) = ' ';
        end
    else                                                                    % Remove anything after char 72, if lines are longer
        if      ismatrix(str),	str(:,w+1:end) = '';
        elseif  ndims(str)==3,  str(:,w+1:end,:) = '';
        end
	end
end





