function str = TABLED1(freq,tableID,xaxis,yaxis,fmin,fmax)

% Set input variable defaults
    if ~exist('tableID','var')	| isempty(tableID)      tableID = 2;            end
 	if ~exist('fmin','var')     | isempty(fmin)         fmin    = 0.0;          end
  	if ~exist('fmax','var')     | isempty(fmax)      	fmax    = 1.e4;         end
    
    if ~exist('xaxis','var') | isempty(xaxis) | regexpi(xaxis,'lin') 	xaxis   = 'LINEAR  ';   end
	if ~exist('yaxis','var') | isempty(yaxis) | regexpi(yaxis,'lin')   	yaxis   = 'LINEAR  ';   end
    if isempty(regexpi(xaxis,'lin')) | regexpi(xaxis,'log')             xaxis   = 'LOG     ';   end
    if isempty(regexpi(yaxis,'lin')) | regexpi(yaxis,'log')             yaxis   = 'LOG     ';   end

% Make freq. and TB column vectors
    freq    = cvec(freq,false);
    nfreq   = length(freq); 
    tableID = cvec(tableID,false);

% Assign appropriate table IDs if necessary
    if numel(tableID)==1
        tableID = [tableID:tableID+nfreq-1]';
    elseif size(tableID)~=size(freq)
        tableID = [2:nfreq+1]';
    end
    
% Generates n by m whitespace char array
    emptyf = @(n,m) repmat(' ',n,m);
    
% Generate the table heading lines
 	L1 = repmat(['TABLED1 ' emptyf(1,8) xaxis yaxis emptyf(1,40)], nfreq, 1);
	L1(:,9:16) = int2field(tableID);

% Second and third line of the tables
	FR = char(num2field(vertcat(fmin,freq,fmax),8));
    F2 = repmat(char(num2field(fmin-0.1)),nfreq,1);
    F3 = repmat(char(num2field(0.)),nfreq,1);
    F7 = repmat(char(num2field(1.)),nfreq,1);
    F32 = repmat(char(num2field(fmax+1.)),nfreq,1);
    F34 = repmat('ENDT    ',nfreq,1);
    
    L2 = [emptyf(nfreq,8) F2 F3 FR(1:nfreq,:) F3 FR(2:nfreq+1,:) F7 FR(3:nfreq+2,:) F3];
    L3 = [emptyf(nfreq,8) F32 F3 F34 emptyf(nfreq,40)];

% Assembly of the tables
	str(1:3:3*nfreq-2,:) = L1;
	str(2:3:3*nfreq-1,:) = L2;
    str(3:3:3*nfreq,:)	 = L3;
    
end


