function [M,dof] = ASM2(M1,M2,dof1,dof2)

%{
    > Assembles M1 and M2 square matrices 
    > Returnes the combined matrix and a combined dof-set
%}

% Sorted column vector of dofs
    dof = unique([cvec(dof1);cvec(dof2)]);
    
% Get Mi to M indices / works with nonsorted Mi
    [~,~,m1m] = intersect(dof1,dof,'stable');
    [~,~,m2m] = intersect(dof2,dof,'stable');

% Expand Mi to M-size and sum
    try 
        M0  = zeros(length(dof));
        M   = zeros(length(dof));
        M0(m1m,m1m) = M1;
        M (m2m,m2m) = M2;  
        M   = M+M0;
    catch
        disp('> M1,M2 must be square matrices <')
    end

% Check if array is a vector and make it column
	function arr = cvec(arr)      
        if any(size(arr)==1)
            if size(arr,1)==1   arr=arr';    end
        else
            disp('> dof must be a cell or a vector <')
            arr = [];
        end       
    end

end

