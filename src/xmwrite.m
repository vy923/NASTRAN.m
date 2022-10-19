%{
-------------------------------------
    Vladimir V. Yotov

    Created: 11.06.2017
    Version: 28.10.2021
-------------------------------------

GENERAL
	Writes matrices for INPUTT4
	Always written as double precision: type = 2/4

	BIGMAT = TRUE format is NOT implemented
	Option to append file not implemented 

	From xmread: [ncols nrows form type name fmt] = nameData{:}
%}

%% 
function xmwrite(M,path,form,name,fmt)

% Capital 'W' works as buffered IO: faster small fprintf calls
fid = fopen(path,'W'); 

% Loop over all matrices to be written
    for j=1:numel(name)
      
    % In case xmwrite is called with only one matrix
        if iscell(M)	Mj = M{j};                
        else        	Mj = M;
    	end
        
	% Matrix type and Fortran to Matlab format conversion
        type(j)   = 2*(~isreal(Mj)+1);
        fmtDat    = sscanf(fmt(j),'%*d%*2c%d%*c%f');                          % "1P,3E23.16" -> [3,23.16]'
        fmtML     = ['\n' repmat(sprintf('%%%gE',fmtDat(2)),1,fmtDat(1))];    % [3,23.16]'   -> '\n%23.16E%23.16E%23.16E'
    
    % Logical masks of nonzero values and columns for writing
        maskVal = (Mj~=0);
        maskCol = any(maskVal,1);
  
    % Print matrix headers
        fprintf(fid,[repmat('%8d',1,4) '%-8s%-s'],flip(size(Mj)),form(j),type(j),name(j),fmt(j)); 
        
 	% Loop over the nonzero matrix columns
        for i=find(maskCol)      
        % Indices of first and last nnz values of a column
            tempVec = maskVal(:,i);
            minInd(i,1)	= find(tempVec,1,'first');
            maxInd(i,1)	= find(tempVec,1,'last');
            
         % Print matrix data line, then all the values of column i
            if type(j)>=3 
                prtVec = Mj(minInd(i):maxInd(i),i);
                fprintf(fid,'\n%8d%8d%8d',i,minInd(i),2*(maxInd(i)-minInd(i)+1)); % --- UNCOMMENT ---
                fprintf(fid,fmtML,[real(prtVec),imag(prtVec)]');                    
            else
            	fprintf(fid,'\n%8d%8d%8d',i,minInd(i),maxInd(i)-minInd(i)+1); % --- UNCOMMENT ---
                fprintf(fid,fmtML,Mj(minInd(i):maxInd(i),i));
            end 
        end
        
   	% Print matrix trailers 
        fprintf(fid,'\n%8d%8d%8d\n%23.16E\n',[1+size(Mj,2) 1 1 1]); % --- UNCOMMENT ---
        
    end % Loop over all matrices 
    
fclose(fid);
end

