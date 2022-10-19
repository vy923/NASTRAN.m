function cleandir(loc,cards)

%{
    > Deletes all files containing 'cards' in the required directory
    > If no 'loc' input exists, current directory is cleaned
    > Does NOT delete old Nastran files, such as f06.1
%}

% Input defaults
	if ~exist('loc','var') | isempty(loc)         
        loc = pwd;        
    end  
	if ~exist('cards','var') | isempty(cards)     
        cards = {'.DBALL','.IFPDAT','.MASTER','.log','.op2','.f04','.asm'};     % op4 usually contains matrices / '.asm'
    end

% File names in the directory
	files = struct2cell(dir(loc));
	files = files(1,:);

% Filter out the files not matching the input cards
	maskf = contains(files,cards);
	files = cellstr(files(maskf));
    
% Delete the mathching files
	for fi = 1:length(files)
    	try     delete([loc '\' files{fi}]);
     	catch   disp('WARNING')
        end
    end
end

