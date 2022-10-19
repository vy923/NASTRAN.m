function str = int2field(vec,padfield,whitespace,N)

%{
    Converts a numeric vector to char array of width N
    Does not catch improper numbers with more than N digits
%}

% Input variable defaults 
	if ~exist('padfield','var') | isempty(padfield)         padfield    = true;  	end  
	if ~exist('whitespace','var') | isempty(whitespace)     whitespace  = false;   	end
    if ~exist('N','var') | isempty(N)                       N = 8;                  end
    
% Convert to int first for faster num2str
    if N>=9     str = num2str(int64(vec),'%-u'); 
    else       	str = num2str(int32(vec),'%-u'); 
    end

% Add whitespace to the left
	if whitespace   
      	str = [repmat(' ',size(vec,1),1) str];
    end
    
% Expand to field width of N, much faster than 'pad' function
	if padfield 
        tmp = repmat(' ',size(str,1),N);
        tmp(:,1:size(str,2)) = str;
     	str = tmp; 
    end
    
end

