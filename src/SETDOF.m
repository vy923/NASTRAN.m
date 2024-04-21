function [dof] = SETDOF(dof,type,scale)

%{
    > Makes a sorted, uniquely indexed degree of freedom list
    
    > Does not catch horizontal cell arrangements / cell spoint
    > Does not catch some noninteger IDs 
    > Cell input must be in format {[grid],[dofs]}
    > String/numerical dof numbers are parsed
%}

% Set default type / scaling
    if ~exist('type','var')  || isempty(type)     type    = 'grid';   	end
    if ~exist('scale','var') || isempty(scale)    scale   = 10;         end   

% Convert the input to dof name vector 
    if  strcmpi(type,'grid')
        
        if iscell(dof)
            for i = 1:size(dof,1)          
                if isempty(dof{i,2})
                    dof{i,2} = (1:6)';
                else
                    dof{i,2} = unique(cell2mat(textscan(num2str(dof{i,2}),'%1f')));
                    dof{i,2} = setdiff(dof{i,2},[0 7:9]);
                end
               	dof{i,1} = setID(cvec(dof{i,1}),dof{i,2},scale);            
            end
        	dof = unique(cell2mat(dof(:,1)));       
        else
            dof = setID(cvec(dof),1:6,scale);
        end  
        
    elseif  strcmpi(type,'spoint')   
        
        try     dof = cvec(unique(scale*dof)); 
        catch   dof = cvec(unique(scale*cell2mat(dof)));
        end 
        
    else
        
        disp('> Type must be GRID or SPOINT <')
        dof = [];       
        
    end 
    
% Check if array is a vector and make it column
	function arr = cvec(arr)      
        if any(size(arr)==1)
            if size(arr,1)==1   arr=arr.';    end
        else
            disp('> dof must be a cell or a vector <')
            arr = [];
        end       
    end

% Takes grid ID and 123456 dofID and creates the ID vector
    function dof = setID(dof,dofID,scale)       
        if ~exist('dofID','var') | isempty(dofID)     
            dofID = 1:6;   	
        end      
    	n   = length(dofID);        
        dof = repmat(unique(dof),1,n);         
     	dof = scale*reshape(dof',n*size(dof,1),1);
        dof = dof + repmat(cvec(dofID),size(dof,1)/n,1);          
    end

end






