function [vecn,dofn,ind,mask,chk] = SELASM(dof,dofSE,lvec)

%{
    > Expands a dof set and load vector to another one
    > Outputs logical mask so dofSE(mask) = intersect(dof,dofSE,'stable')
    > Outputs an index matrix with [dof dofSE] indices of the intersection
    > Parameter chk == 1 if dof is entirely part of dofSE
    > Assumes dof vectors with SETDOF format. 
    > NOTE: No conjugate transpose of complex vector inputs is done in the function (check?)
%}

% Set defaults / check if lvec and dof are same length
    if      ~exist('dofSE','var')| isempty(dofSE)       dofSE = [];       end
  	if      ~exist('lvec','var') | isempty(lvec)        lvec  = [];  
    elseif  any(size(dof)~=size(lvec))                  return
    end

% Find intersection indices
    [~,ind,indSE] = intersect(dof,dofSE,'stable');
	chk = (length(indSE)==length(dof));

% Initialise dofSE-sized lvec, dof
    vecn = zeros(size(dofSE));
    mask = zeros(size(dofSE));

% Logical index mask on dof / index matrix / remaining dof set
	mask(indSE) = dof(ind);
    mask = mask~=0;
    dofn = dof(ind);
    
    try     vecn(indSE) = lvec(ind); 
    catch   vecn = [];   
    end

% Combine indcices 
    ind	= [ind indSE];

end

%{
% Example of assembling two complex loads onto dofSE

    lvec1 = 1./(1-2*rand(4,1)) + sqrt(-1)./(1-2*rand(4,1));
    lvec2 = 1./(1-2*rand(11,1)) + sqrt(-1)./(1-2*rand(11,1));

    dof1    = SETDOF({[1 6 11 16],1});
    dofSE	= SETDOF({[11 1 20 5 40 50],1});
    dof2    = SETDOF({[11 1 20 30 40 16 4 7 9 6 55],1});

    [vecn1,dofn1] = SELASM(dof1,dofSE,lvec1)
    [vecn2,dofn2] = SELASM(dof2,dofSE,lvec2)

    lvec_tot = vecn1+vecn2
%}
