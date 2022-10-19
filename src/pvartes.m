function [A,qdof] = pvartes(adof,qdof,seID,pathout)

%{
    > Works ONLY with 123456 
%}

% To avoid a lot of slow num2str
    charnum(1:7) = ['1','2','3','4','5','6','0'];
    
% Complies with EXTSEOUT format [9sssnnnn]
    qdof = 9e7 + 10e3*seID + qdof;

% Take the format of the two header lines 
    A{1} = 'DMIG    matnameI       0       9       2       0                       1';
    A{2} = 'DMIG*           matnameI               1               0                *';
    line = '*               pointnum               0    1.00000D+00                 *';

% Name the matrices SE[seID] and update header
    namestr = pad(['SE',sprintf('%d',seID)],8,'left');
    A{1}    = strrep(A{1},'matnameI',namestr);
    A{2}    = strrep(A{2},'matnameI',namestr);
    
% Write a-set to A
    for i = 1:length(adof)
        for j = 1:6
            lnum        = 2+6*(i-1)+j; 
            temp        = line; 
            temp(40)    = charnum(j);
            temp(17:24) = sprintf('%-8d',adof(i));
            A{lnum}     = temp;
        end
    end
    
% Write q-set to A and remove the last line continuation sign
    for i = 1:length(qdof)
        temp        = line; 
        temp(17:24) = sprintf('%-8d',qdof(i));
        A{lnum+i}   = temp;
    end
        A{lnum+i}(length(line))=[];
        
% Write A to the output file
    textout(A,pathout);

end

