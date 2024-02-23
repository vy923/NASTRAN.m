function [R,dat] = resread(A,card,quickSearch,emptyval)

%{
    > Scans an f06 file and returns the type of output requested on 'card'
    > A can be a cell vector w/o EOL chars, as returned by textin(X,false)
    > If the above is true, A = string(A): improved textscan/contains speed
    > quickSearch improves file scanning speed, especially for large files

    > Results are returned as R{subcase}(frequency,component,pointID)
    > Eigenvectors are returned as R{subcase}(pointID,value,modeNumber)
    > dat(subcase,:) = {freq;nodes}. For eigenvectors dat{subc} = modeNumber
    
    > Block reading offset 'off(1:3)' = [card to blockStart, flagend to 
    blockEnd, pageval to blockStart]
    > blockEnd may have a higher offset w.r.t flagEnd due to some Nastran
    text output just after the numeric data, but textscan ignores this
    > 'pageval' contains PointID or Frequency, and Cycles and Mode Number
    if real eigenvectors are read

    > Assumed behaviour as in Nastran v2012.1: Only one output card of a 
    given type per subcase is allowed, if card is repeated, only the last 
    one is output. Also, one SORTi value per subcase is allowed, SORT2 
    takes prevalence if both are present. 

    > flagEnd = 'MSC.NASTRAN' for 2012.1. Currently set for v2017/v2020

CHANGES
    - MPC/SPC card ' ' offset adjusted 25 -> 27 for v2020
    - subcase calc: blockEnd+off(2)+2 -> blockStart-off(3)-1 for v2020

PENDING
    - (1.0) fix chk=2 / static solutions not working, spec. pageval calc
    - DOESN'T WORK FOR PHASE...
%} 

    if regexpi(card,'disp')     chk = 1;	card = [repmat(' ',1,39) 'C O M P L E X   D I S P L A C E M E N T   V E C T O R'];                            	end
    if regexpi(card,'acc')      chk = 1;    card = [repmat(' ',1,39) 'C O M P L E X   A C C E L E R A T I O N   V E C T O R'];                           	end
    if regexpi(card,'oload')	chk = 1;    card = [repmat(' ',1,47) 'C O M P L E X   L O A D   V E C T O R'];                                           	end
    if regexpi(card,'spcC')     chk = 1;    card = [repmat(' ',1,27) 'C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T'];	end
    if regexpi(card,'mpcC')     chk = 1;    card = [repmat(' ',1,27) 'C O M P L E X   F O R C E S   O F   M U L T I P O I N T   C O N S T R A I N T'];      end
    if regexpi(card,'vel')      chk = 1;    card = [repmat(' ',1,43) 'C O M P L E X   V E L O C I T Y   V E C T O R'];                                    	end
    if regexpi(card,'vec')      chk = 2;    card = 'R E A L   E I G E N V E C T O R   N O .';                                                               end
	if regexpi(card,'spcR')     chk = 2;    card = 'F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T';                                     end
	if regexpi(card,'mpcR')     chk = 2;    card = 'F O R C E S   O F   M U L T I P O I N T   C O N S T R A I N T';                                         end
    if regexpi(card,'txdd')     chk = 2;    card = [repmat(' ',1,44) 'A C C E L E R A T I O N    V E C T O R'];                                             end

% Default replacement of empty values / quickSearch request    
    if ~exist('emptyVal','var')     	emptyval    = NaN;          end
    if ~exist('quickSearch','var')    	quickSearch = true;         end  
    if quickSearch~=1 | chk==2      	quickSearch = false;        end         
    
% Convert A to string to accelerate execution of some funcitons
    if iscell(A)    A = string(A);
    end

% Exit function if card input is inappropriate
    if ~exist('chk')
        disp('WARNING: Unsupported output request specified on card')
        R = double.empty;
        return
    end
    
% Formats and offsets for different cards
    if chk == 1
      	off = [4,1,5];
        fmx = ['%*s%*s%f'];
        fmt	= ['%*f%f%*s' repmat('%f%f%f%f%f%f%*[^\n]\n',1,2)]; 
    elseif chk == 2
     	off = [4,1,4];                                                      % MOD 27.10.21 [3,1,3] -> [3,1,4] / 12.02.24 -> [4,1,4]
        fmx = ['%*s%*s%f' repmat('%*s',1,sum(char(card)~=' ')) '%f'];
     	fmt	= ['%f%*s' '%f%f%f%f%f%f'];
    end

% End of block/sort1/sort2 flags
    if quickSearch  
        flagEnd = string(['1' repmat(' ',1,70)]);  
    else
        flagEnd	= "MSC Nastran";                                            % MSC.NASTRAN for 2012.1
    end
    flagS1  = "FREQUENCY";
    flagS2  = "POINT-ID";
    flagMP  = "MAGNITUDE";                                                  % To identify (REAL/IMAG) or (MAGNITUDE/PHASE) output
    card    = string(card);                                                 % Char search is faster, but quickSearch requires a string

% Find the initial and final line of each block
    if quickSearch
        tmp         = find(startsWith(A,[card,flagEnd]));                   % Was char cell: {card,flagEnd}. Slower.
     	blockStart  = tmp(startsWith(A(tmp),card)) + off(1);
    	blockEnd    = tmp(startsWith(A(tmp),flagEnd)) - off(2);
    else
        tmp         = find(contains(A,[flagEnd strtrim(card)]));
        blockStart  = tmp(contains(A(tmp),strtrim(card))) + off(1);     
        blockEnd    = tmp(contains(A(tmp),flagEnd)) - off(2);
    end

% Terminate if no result blocks are found
    if isempty(blockStart)
       disp('NOTE: No results match the requested card. Try quickSearch = false')
       R    = cell.empty;
       dat  = double.empty;
       return
    end

% Final line of each block; matrix is logical, so 1 byte per entry 
	[~,ind]  = max(blockEnd>blockStart');
	blockEnd = blockEnd(ind)
    
% Check if block coordinates are correct
    if length(blockEnd)~=length(blockStart)
        disp('WARNING: Inappropriate block termination flag')
       	R = cell.empty;
        return
    end

% Sort value, subcase number, pageval = Point/Freq/eigenfr.+modeNumber
    sortval = contains(A(blockStart-off(1)-1,:),flagS2) + 1;
    pageval = cell2mat(textscan(join(A(blockStart-off(3),:)),fmx));
	subcase = cell2mat(textscan(join(A(blockStart-off(3)-1,:)),'%*u%*s%f'));    % MOD 17.09.21
    complx  = contains(A(blockStart-off(1)+1),flagMP);                          % Real/Imag = false, Mag/Phase = true

	if chk==2                                                               % MOD 28.10.21: temporary fix of (1.0)
        pageval = subcase;  
        tmpPar = 2; 
    end
    
% Read all identified blocks 
    for i = 1:length(blockStart)
        As = join(A(blockStart(i):blockEnd(i)),newline);
        try     S{i} = cell2mat(textscan(As,fmt,'CollectOutput',0));
        catch   S{i} = cell2mat(textscan(As,fmt,'CollectOutput',1,'EmptyValue',emptyval));
        end
    end

% Function converting magnitude/phase to complex
    mp2complex = @(mag,phase) mag.*exp(sqrt(-1)*deg2rad(phase));

% Assemble the blocks
    try
    if chk==1
        subcs = unique(subcase);
        maskS = sortval==1;    
        
    % Loop through all the subcases
        for i = 1:length(subcs)  
            maskS1 = (subcase==subcs(i) & maskS);
            maskS2 = (subcase==subcs(i) & ~maskS);

     	% Different assembly process depending on SORTi value
            if any(maskS1)
                tmpID = unique(pageval(maskS1));
                for j = 1:length(tmpID)
                    Q(:,:,j) = cell2mat(S(maskS1 & pageval==tmpID(j)));
                end  
              	node = Q(:,1,1);
              	freq = tmpID;
               	Q = permute(Q,[3 2 1]); 
            else
                tmpID = unique(pageval(maskS2));
                for j = 1:length(tmpID)
                    Q(:,:,j) = cell2mat(S(maskS2 & pageval==tmpID(j))');
                end
              	node = tmpID;
             	freq = Q(:,1,1);
            end
            
      	% Collect the assembled result for i-th subcase
            if any(complx(maskS1|maskS2))                                   % All values corresponding to the subcase
                R{i} = mp2complex(Q(:,2:7,:),Q(:,8:13,:));
            else
                R{i} = Q(:,2:7,:) + sqrt(-1)*Q(:,8:13,:);
            end
                dat(:,i) = {freq;node}; 
                    clearvars Q freq node  
        end     
    
    elseif chk==2                                                           % MODDED - 22/1/19: LINE limit number is insufficiently high and responses are split in blocks
        subcs = unique(subcase);  
        for i = 1:length(subcs)    
            [tmpID,idxTMP] = unique(pageval(subcase==subcs(i),2/tmpPar));   % MOD A: Added idxTMP to store indices of unique pageval values
            for j = 1:length(tmpID)  
              	Q(:,:,j) = cell2mat(S(subcase==subcs(i) & pageval(:,2/tmpPar)==tmpID(j)).');   % MOD B: Added transpose (.'). Works for split and full result blocks. 
            end
         	R{i} = Q(:,2:7,:);
            dat{i,1} = [tmpID pageval(idxTMP,1)];                           % MOD C: added pageval(idxTMP,1), contains frequencies when reading 'vec'. CHECK for other cases!
            dat{i,2} = Q(:,1);                                              % MOD D: return node numbers
                clearvars Q 
        end
        
    end %Subcase loop
        
    catch
        disp('WARNING: Results assembly failed, R returned in cell format');
        R = S;
    end

end   



