function str = num2field(arr,maxd)

%{
    > Writes a matrix arr onto a numel(arr)-by-maxd char vector
    > Maximum possible precision is retained, depending on maxd
    > Exponent sign is removed, as used by Patran/Nastran (e.g. 1.7-5)
    > Rounding is done by sprintf, no loss of precision due to cut-off
    > UPDATE [22/3/2018]: Fixed a bug in maskS1 case
    > UPDATE [30/06/22]: check if input is not a vector and expand

    > Optional / incomplete
        Var inputs
        Reshape if input is not a vector
        Calling with int16/32/64
        Engineering notation output
%}

% Initialisations / incomplete
 	if ~exist('maxd','var') || maxd<6       maxd = 8;     end               % Default printed field width
    if ~isvector(arr), arr=arr(:);                        end
    
% Exponents of the numbers in engineering notation
 	expn = floor(log10(abs(arr)));
    
% Masks for NaN, zero and small exponents
    maskN   = isnan(expn);                                                  % NaN
    maskF   = (expn>=-2 & expn<=maxd-3 & ~maskN);                           % Small exponent
  	maskFE  = (maskF & expn>=0);                                            % Numbers potentially needing correction

% Correction for rounding in sprintf resulting in maxd+1 field outputs
    if any(maskFE) 
        maskT = (abs(arr(maskFE)).*10.^(maxd-3-expn(maskFE)) >= 10^(maxd-2)-0.5);     
        if any(maskT)
            indTmp = find(maskFE);
            indTmp = indTmp(maskT);
            arr(indTmp)  = arr(indTmp) + 0.5/10^(maxd-2);                   % Same effect as round(..., 'significant');
            expn(indTmp) = expn(indTmp) + 1;
            
            q = expn(indTmp)>maxd-3;                                        % In case after correction maskF overlaps with mask S1
            if any(q)
                maskF(indTmp(q)) = false;
            end
        end
    end

% Remaining masks
  	mask0   = (arr==0);                                                     % Zero
  	maskS1	= ~(maskF | mask0 | maskN | abs(expn)>=10);                     % Scientific, 1d exponent
  	maskS2	= ~(maskF | mask0 | maskN | abs(expn)<10);                      % Scientific, 2d exponent
    maskFN  = (maskF & expn<=0);                                            % Non-positive exponents of small exp. numbers
   
% Initialisation of the output
  	str = repmat(' ',numel(arr),maxd);                                      % Optional

% Print most values onto temporary character vectors
    if any(mask0)
        str(mask0,:) = repmat(pad(' 0.',maxd,'right'),nnz(mask0),1);
    end

    if any(maskFN)
        fmtFN = sprintf('%% %.0f.%.0ff',maxd,maxd-3);                       % Produces '% 8.5f'
        strFN = sprintf(fmtFN,arr(maskFN)); 
        str(maskFN,:) = reshape(strFN,maxd,numel(strFN)/maxd)';
    end

    if any(maskS1)
        fmtS1 = sprintf('%% %.0f.%.0fe',maxd,maxd-5);
        strS1 = sprintf(fmtS1,arr(maskS1));
        strS1 = strrep(strS1,'e+0','+');
        strS1 = strrep(strS1,'e-0','-');
        strS1 = strrep(strS1,'0e+1','+1');                                  % UPDATE [22/3/2018]: fixes the case  9.99999e+9 -> '1.000e+10' and strrep does not work
        str(maskS1,:) = reshape(strS1,maxd,numel(strS1)/maxd)';
    end

    if any(maskS2)
        fmtS2 = sprintf('%% %.0f.%.0fe',maxd,maxd-6);                        
        strS2 = sprintf(fmtS2,arr(maskS2)); 
        strS2 = erase(strS2,'e');
        str(maskS2,:) = reshape(strS2,maxd,numel(strS2)/maxd)';
    end

% Special treatment for small exponent values > 0, maximising sig. digits 
    for i = 1:maxd-3
        maskFx = (expn==i);        
        if any(maskFx)   
            if i==maxd-3    
                fmtFx = sprintf('%% %.0f.%.0ff.',maxd-1,maxd-3-i);
            else
                fmtFx = sprintf('%% %.0f.%.0ff',maxd,maxd-3-i);
            end   

            strFx = sprintf(fmtFx,arr(maskFx));
            strFx = reshape(strFx,maxd,numel(strFx)/maxd)';
            str(maskFx,:) = strFx;
        end     
    end

end

%{
Test for some special cases: 
for i=-20:20
    disp(num2field(-9.9999999999*10^i,8))
end
%}



