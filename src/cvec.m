%{
-------------------------------------
    Vladimir V. Yotov
    Te Pūnaha Ātea Space Institute
    University of Auckland

    Version: 16.12.2021
-------------------------------------

NOTES
	Check if array is a vector and make it column
	clr = true makes arr empty on return
	opt = false returns error if arr is not vector
    ***Relatively high call overhead***
%}

function arr = cvec(arr,opt,clr)

% Initialisations
  	if ~exist('opt','var')      opt = 1;       	end
  	if ~exist('clr','var')      clr = 0;       	end  
 
% Main function
    if ismatrix(arr) && any(size(arr)==1) 
        if size(arr,1)==1   arr = arr.';
        end
    elseif opt
      	disp('> Input is not a vector <')
        if clr 	arr = [];
        end
    else
        error('> Input is not a vector <')
    end  

end

