function str = DLOAD(dof,vec,vectype,option,startID,dloadID,tableID,scale)

%{
    > Input 'dof' is assumed in SETDOF format for non-pressure loads.
    > Input 'dof' is assumed element ID for the case of pressure loads.
    > Input 'dof' must be a vector, constant for all columns of 'vec',
    thus multiple loads allowed, but only to one set of DOFs

    > SID/EIDL/ELEM ID sets can have overlapping IDs.
    > SID values must be unique and different from dloadID: this is checked
    > It is assumed that each load vector contained in columns of 'vec'
    will reference a separate TABLED. One tableID per frequency is defined. 
    
    > Option 0 assigns a load card to each loaded DOF at each freq. and     
  	unit scaling in DLOAD
    > Option 1 defines a unit load for all loaded DOFs and scales them to
 	the appropriate magnitude in DLOAD

    > NOTE: If load cards FORCE/MOMENT/PRESSURE/SLOAD have repeated IDs, 
    their values are first summed over the repeated ID, after which the
    total value is used by RLOAD cards. 
    > NOTE: If several RLOADs reference the same load card, the overall 
    load due to calling that card multiple times is summed and scaled 
    according to the applied phases in RLOAD and applied scaling in the 
    referenced TABLEDs and DLOAD. 

    + MOD/FIX 18.04.21 to work with input load of less than 8 DOFs
%}

% Set input variable defaults
    if ~exist('vectype','var') | isempty(vectype)       vectype	= 'force'; 	end     % Force/moment/sload or pressure/pload
    if ~exist('option','var')  | isempty(option)        option	= true;    	end     % Optimised output format if true
 	if ~exist('dloadID','var') | isempty(dloadID)       dloadID = 1;       	end   	% DLOAD card must differ from RLOAD IDs
  	if ~exist('startID','var') | startID <= 0           startID = [];      	end   	% If ~empty, SID/EID = startID:startID+n-1 
    if ~exist('tableID','var') | numel(tableID)~=1      tableID = 2;       	end    	% Dynamic load table ID
 	if ~exist('scale','var')   | isempty(scale)         scale   = 10;      	end   	% As dof should be in SETDOF format / irrelevant for pressure load

% Type of input/required output. 
    if      regexpi(['force','moment','sload'],vectype) 	vectype = 1;    % if(double.empty) is identical to if 0    
    elseif  regexpi(['pressure','pload'],vectype)           vectype = 2;       
    end                          
   
% Exit function if vectype input is inappropriate    
    if ~ismember(vectype,1:2)
        disp('WARNING: No output requests match vectype')
        str = char.empty;
        return
    end
    
% Expand the inputs to n-by-q size
    if length(vec)~=numel(vec)
      	q = size(vec,2);
        try     
            dof = repmat(cvec(dof,false),1,q);                              % cvec catches empty/improper shape arrays                           
        catch
            disp('WARNING: dof must have size n or n x q')
          	str = char.empty;
         	return 
        end  
    else
        q = 1; 
    end

% Assigns unique table IDs for each separate load (column) in 'vec'    
    try     tableID = tableID:tableID+q-1;
    catch   tableID = 2:q+1;
    end
    
% Remove zero load values and reshape inputs to column vectors
    mask = vec~=0;
    vec  = cvec(vec(mask));                                                 % Always returns column vectors                                            
    dof  = cvec(dof(mask)); 
    n    = length(vec); 
    
% Adjust the input variables based on vectype/q
	if vectype == 2
        scale = 1;                                                          % Stops division of the element IDs by scale 
        disp('NOTE: scale = 1 assigned')
    end                
    if isempty(startID) && q>1                                              % Prevents generation of non-unique SID values
        startID = 2;
        disp('NOTE: startID = 2 assigned')
    end                 

% Variables containing the loaded nodes/elements
    if option 
      	[uniD,~,indD] = unique(dof);                                        % EID references ELEM(indD)
        ELEM = floor(uniD/scale);              
    else
    	ELEM = floor(dof/scale);     
    end

% Ensure that RLOAD2 card IDs are unique and do not clash with dloadID
    if ~isempty(startID) || q>1     
        SID = [startID:startID+n-1]';   
    else
        SID = dof;
    end 
    if any(SID==dloadID)
        disp('WARNING: RLOAD IDs overlap with the requested dloadID')
        str = string.empty;
        return
    end

% Check if ELEM/SID numbers fit on 8 fields and convert to char
    if any(SID>=10e7) || any(ELEM>=10e7)
        disp('WARNING: Loaded node/element IDs are larger than 8 digits')
        str = string.empty;
        return
    else
        ELEM = int2field(ELEM);
        SID  = int2field(SID);
    end
    
% Utility functions for generating empty fields, ones and zeros
    emptyf = @(n,m) repmat(' ',n,m);
    nulstr = @(n,m) repmat(num2field(0),n,m);                               % Obsolete: char(arr2field(0))
    onestr = @(n)   repmat([num2field(1,7) ' '],n,1);                       % Obsolete: char(arr2field(1,7))
    
% Transform the load vector to magnitude and phase char array
vec
    MAGN  =  num2field(abs(vec),8);                                         % Obsolete: char(arr2field(abs(vec),7));
    PHASE =  num2field(rad2deg(angle(vec)),8);                              % Obsolete: char(arr2field(rad2deg(angle(vec)),7));    
    
% Table ID cards, num2str calls somewhat minimised
	TB = string(int2field(tableID',1,1))';
 	TB = repmat(TB,n/q,1);
	TB = char(reshape(TB,n,1));
 
% Load card number/called ID and scaling factors
 	if option
        disp('NOTE: Load IDs EIDL = ELEM when option = true')
%        EIDL    = SID(1:length(ELEM),:);                                    % Was EIDL=ELEM, giving repeated EIDL for all dofs of a node      
        EIDL    = SID(1:size(ELEM,1),:);                                    % MOD 18.04.21 / length(ELEM)->size(ELEM,1), otherwise does not work for lvec size < 8
        EIDRL   = EIDL(indD,:);                                             % Was ELEM(indD,:);
        DLSCL   = MAGN;                                                     % All magnitudes defined in DLOAD
    else
        EIDL    = SID;
        EIDRL	= SID;
        DLSCL	= onestr(n);
    end 
    
% Generate RLOAD2 cards
    RLOAD = padfield([repmat('RLOAD2  ',n,1) SID EIDRL emptyf(size(EIDRL,1),8) PHASE TB]);
   
        clearvars mask vec PHASE TB indD
    
% Generate the force/moment/sload cards
    if vectype == 1     
        
    % When fast option is used, no repeated unit loads are allowed
      	if option   temp = rem(uniD,scale);
        else        temp = rem(dof,scale);
        end   
        
    % Split the load into force/moment/sload
        maskF = (temp>=1 & temp<=3);
        maskM = (temp>=4 & temp<=6);
        maskQ = (temp==0);    
        
    % Column of Nastran cards
      	CARD(maskF,:) = repmat('FORCE   ',sum(maskF),1);                    % Alternatively find/logical(mask)
        CARD(maskM,:) = repmat('MOMENT  ',sum(maskM),1);
        CARD(maskQ,:) = repmat('SLOAD   ',sum(maskQ),1); 
        
    % Fields 4,5 for the above cards
        F4 = emptyf(length(temp),8);                                        % CID field of force/moment or unit SLOAD
        F5 = emptyf(length(temp),8);                                        % Total scaling of force/moment, empty for SLOAD
        
     	if option
            F4(maskQ,:)     = onestr(sum(maskQ));      
            F5(~maskQ,:)    = onestr(sum(~maskQ));
        else
            F4(maskQ,:)     = MAGN(maskQ,:);
        	F5(~maskQ,:)    = MAGN(~maskQ,:); 
        end
   
    % Temporary variables
        remtemp	= rem(temp,3);                                              % Used to select N1,N2,N3 dof of force/moment
        remtemp(maskQ) = NaN; 
        onetemp	= onestr(length(remtemp));                                  % Vector of ones for the nonzero N1,N2,N3 values 

 	% Fields 6 to 8 of the above cards    
        FN = emptyf(length(temp),24);                                       % Force/moment vector components, empty for SLOAD
        FN(~maskQ,:) = nulstr(sum(~maskQ),3);                               % The rest of FN is initialised to 0.0   
    
  	% Nonzero parts of FN
     	FN(remtemp==1,1:8)   = onetemp(remtemp==1,:);                       % Make nonempty N1 fields = 1.0
      	FN(remtemp==2,9:16)  = onetemp(remtemp==2,:);
      	FN(remtemp==0,17:24) = onetemp(remtemp==0,:);  
        
        CARDS = padfield([CARD EIDL ELEM F4 F5 FN],72);   
    end
    
% Generate the pload4 cards if requested    
    if vectype == 2
        CARD = repmat('PLOAD4  ',size(EIDL,1),1);
        if option   
            CARDS = padfield([CARD EIDL ELEM onestr(length(uniD))],72); 
        else
            CARDS = padfield([CARD EIDL ELEM MAGN],72);
        end     
    end

        clearvars *temp* mask* CARD EIDL ELEM MAGN 
        clearvars -regexp ^F\S{1}$

% Assemble the DLOAD part of the output file: independent of vectype        // CAN be optimised with padfield     
    nrows   = ceil((n-3)/4+1);
	nblanks	= rem(4-rem(n+1,4),4);
  	DLOAD   = [int2field(dloadID), onestr(1); DLSCL, SID; emptyf(nblanks,16)]; 
  	DLOAD   = reshape(DLOAD',64,size(DLOAD,1)/4)';
 	DLOAD   = [vertcat('DLOAD   ',emptyf(nrows-1,8)) DLOAD];  
    
        clearvars DLSCL SID
        
    str = vertcat(CARDS,RLOAD,DLOAD);
    
end

%	[RLOAD2 SID  EIDRL emptyf(n,8) PHASE TB];
%	[CARD   EIDL ELEM  emptyf(n,8) LSCL];

% ----------------------------- FUNCTIONS ---------------------------

% Faster version of 'pad', returns char  
    function str = padfield(str,width)                                      % Slightly faster than direct indexing of str
        if ~exist('width','var')    
            width = 72;      
        end         
        if size(str,2) == width
            return
        else
            tmp = repmat(' ',size(str,1),width);
            tmp(:,1:size(str,2)) = str;                                     % Aligns to the left
            str = tmp;
        end
    end

% ------------------------------ EXAMPLES ---------------------------    
%{
    DLOAD(SETDOF([888,999]),cvec([1,1,1,0,0,0,2,2,2,4,0,7]),[],true,[],2,1,[])    
    DLOAD(SETDOF(999),cvec([1,1,1,0,0,0]),[],[],[],2,1,[])  
%}




