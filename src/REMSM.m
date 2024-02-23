function [M,SMset,dof] = REMSM(M,dof,SMset,type,scale)

%{
    > Name originates from nastran S- and M-sets 
    > MPC format should be {[master],[dependent],[component]}
    > MPC{i,:} can have k/n master/slave nodes, k must be 1 or n            //	No verification for is k implemented
    > DOF input scaling should be consistent with 'scale'
    > M should be a vector or square matrix                                 //	Input array dimensions are not fully checked
    
    > if M is a load matrix, it must not be square for 'spc' type
    > if M is a load matrix, function SEEMS TO work fine for 'mpc'

    UPDATE [22.4.21] 
    > Updated REMSM to check if dof-removing while loop gets stuck
    > To avoid infinite loops, another check was added s.t. in case of 
    nested RBE2s only lowest level ('true') slave DOFs are removed at each 
    stage. 
    > Cyclic RBE2s should not work

    NOTE: 
    - Define a multilevel 'RBE2 tree' and test if REMSM works?
%}

% Set some argument defaults
    if ~exist('scale','var')	scale   = 10;           end
    if ~exist('type','var') 	type    = 'SPC';    	end
    
% Set type check
 	if      strcmpi(type,'SPC')   	spc = SMset; mpc = [];
    elseif  strcmpi(type,'MPC')  	mpc = SMset; spc = []; 
    else
        disp('> Type must be SPC or MPC <'); return
    end
    
% Check if M is a vector or square matrix
 	try     dm=1; M=cvec(M,false); 
 	catch	dm=2;   
    end
     
% Check if M is a non-square load matrix                                    UPDATED: [29/3/2018]
    if strcmpi(type,'MPC') && dm==2 && size(M,1)~=size(M,2)
        disp(':: REMSM :: M not square, ASSUMED to be a load matrix of size [dof,nfreq]');
        dm = 1;
    end

% Process SPC first to speed up mpc removal
    if ~isempty(spc)
        spc         = SETDOF(spc,'grid',scale);
        [dof,index]	= setdiff(dof,spc,'stable');
        
        if (dm==1)	
            M = M(index); 
        else
            try     M = M(index,index);
            catch   M = M(index,:);                                         % NOTE: DOES NOT CATCH A SQUARE LOAD MATRIX INPUT
            end
        end
        
        SMset = spc; 
    end
   
% Process RBE2 MPC generating slaveDOF to masterDOF vector
    if ~isempty(mpc)
        for i = 1:size(mpc,1)
            try
                mpc{i,1} = SETDOF(mpc(i,[1 3]),'grid',scale);
                mpc{i,2} = SETDOF(mpc(i,[2 3]),'grid',scale);
                mpc{i,1} = repmat(mpc{i,1}, length(mpc{i,2})/length(mpc{i,1}), 1);
            catch
                disp('> Only one master grid per RBE2 is allowed <')
                mpc = [];
                return
            end
        end
        mpc = unique(cell2mat(mpc(:,[1 2])),'rows');

    % Remove MPC dofs that are not in M's dof set
        mpc     = mpc(all(ismember(mpc,dof)'),:);                           % Both master and dependent are in 'dof'
        SMset   = mpc; 

    % Summation of the dependent onto independent mpc dofs
        while true   

            [~,indMPC,indMm] = intersect(mpc(:,1),dof,'stable');            % mpc set subindex / master to M
            
            mask    = ~ismember(mpc(indMPC,2),mpc(:,1));                    % MOD [22.4.21] Remove master DOFs that are also slave in a higher level RBE2
            indMPC  = indMPC(mask);                                         % hence only 'true' slave DOFs are processed at each stage    
            indMm   = indMm(mask);
            
            [~,~,indMs] = intersect(mpc(indMPC,2),dof,'stable');            % dependent to M
            
            % Output vector/matrix processing
                M(indMm,:) = M(indMm,:) + M(indMs,:);                       % NOTE [20.4.21] Irrelevant when removing RBE2s
                M(indMs,:) = [];               
            if (dm~=1)	
                M(:,indMm) = M(:,indMm) + M(:,indMs);
                M(:,indMs) = [];
            end

            testVar        = size(mpc);                                     % MOD [22.4.21] added to break infinite while loops
            mpc(indMPC,:)  = [];
            dof(indMs)     = [];     
            
            if isempty(mpc)                                                 % End loop
                break    
            elseif all(size(mpc)==testVar)                                  % MOD [22.4.21] added to break infinite while
                error(':: REMSM :: Infinite while loop detected')
            end
                         
        end % WHILE outer
        
    end % IF MPC

end





