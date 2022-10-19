function [vec] = str2val(arr,emptyval)

%{
    > Reads numbers from a vector written on an n x m char array
    > One number per line is expected, possibly in Patran short format
    > By default empty fields are returned as NaN
    
    > Can be called with string/char/cell vectors
%}
  
% Return in case of empty input
    if  isempty(arr)
        vec = arr;  
        return      	
    end

% Set varable defaults
    if ~exist('emptyval','var')     emptyval = NaN;             end
    if ~ischar(arr)                 arr = char(arr);            end
    
% Prepare array for reading
	arr = arr';
	ind = ~all(arr==' ')';
	arr(size(arr,1)+1,:) = newline;
    
% Initialise the output vector
    vec = emptyval*ones(size(ind));
    
% Construct the final float vector
    if any(ind)
        tmp = cell2mat(textscan(arr(:,ind),'%f%f','EmptyValue',0,'CollectOutput',1));
      	tmp = tmp(:,1).*10.^tmp(:,2);
        vec(ind) = tmp;
    end
    
end


%{
% --------------------------- OLD VERSION ---------------------------
% Return in case of empty input
    if  isempty(arr)
        vec = arr;  
        return      	
    end

% Set varable defaults
    if ~exist('emptyval','var')     emptyval = NaN;             end
    if ~isstring(arr)               arr = string(arr);          end

% Treatment depends on the presence of trailing whitespaces 
    if all(strlength(arr)==strlength(arr(1)))
        ind = arr~=repmat(' ',1,strlength(arr(1)));
    else
        ind = arr~="";                                                      % arr = erase(arr," "); % Optional
    end   

% Initialise output to the empty value
	vec = emptyval*ones(size(ind));   

% Construct the final float vector
    if any(ind)
        tmp = cell2mat(textscan(join(arr,newline),'%f%f','EmptyValue',0,'CollectOutput',1));
        tmp = tmp(:,1).*10.^tmp(:,2);
       	vec(ind) = tmp;
    end
%}

