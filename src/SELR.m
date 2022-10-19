function [FT,FQ] = SELR(Lvec,Lset,Gset,Tset,Qset,GOT,GOQ)

%{  
    Inputs lvec, lset contain the dof-load pairs
    Input transoformation matrices GOT/phiR, GOQ/phiL
  	Assumes dof vectors with SETDOF format, also implying sorted
    G-set may or may not include Q-set
%}

% Define O-set
    %Oset = setdiff(Gset,[Qset;Tset]);                                     % Original version
    Oset = Gset(~ismember(Gset,[Qset;Tset]));                               % UPDATE [28/3/18]: Allows repeated entries in the Gset
    
% Split the load to T/O-set: empty dof means vec == 0   
    if any(size(Lset)~=size(Lvec))                                          % n x nfreq Lvec input matrix
        [~,~,indT] = SELASM(Lset,Tset);
        [~,~,indO] = SELASM(Lset,Oset);
        vecT = Lvec(indT(:,1),:);
        vecO = Lvec(indO(:,1),:);
    else                                                                    % n x 1 Lvec input: single load
        [vecT,dofT,indT] = SELASM(Lset,Tset,Lvec);
        [vecO,dofO,indO] = SELASM(Lset,Oset,Lvec);                          % vecO = Lvec(indO(:,1))
    end
        
disp([size(vecO) size(vecT) size(GOT) size(GOQ)])
% Reduce the T/O-set loads
    try
        FT = GOT.'*vecO + vecT;
        FQ = GOQ.'*vecO;
    catch
        disp('> Incorrect set definitions or matrix sizes <')
        return
    end
    
end

%{ 
DOF sets, nastran notation
	Q: Modal reduced/physical
	T: Reduced physical
	A: T+Q / Total reduced
	O: Omitted
	G: A+O / T+O+Q / Total physical   
%}

%{
% Example 
    Lvec    = 1./(1-2*rand(11,1)) + sqrt(-1)./(1-2*rand(11,1));
    Lset	= SETDOF({[11 1 20 30 40 16 4 7 9 6 91],1});
    
    Tset    = SETDOF({110*[1 4 7 9 6],1});
    Gset    = SETDOF({[1 4 7 9 6 30 40 20 16],1});
    Qset    = SETDOF({91 92 99},'spoint');
%}
