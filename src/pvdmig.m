%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       PVDMIG generates a DMIG input with DOF lists for manual CB reassembly
%
%   See also: pvartes
%
%   NOTES
%       [0 9 1 0 NaN NaN 1] -> [0 9 2 0 NaN NaN 1] double prec
%
%   VERSION
%   v1.0 / 16.02.23 / V.Y.
%  ------------------------------------------------------------------------------------------------

function A = pvdmig(adof,qdof,seID)

% Numeric DMIG entries
    tmp = repmat(adof(:)',6,1);
    tmp = [ tmp(:) repmat(cvec(1:6),numel(adof),1).*[1 0] + [0 1];
            qdof(:).*[1 0 0] + [0 0 1] ];

% Header for real 1 col matrix
    A = ['DMIG    ' pad(seID,8,'right') sprintf('%-8u',[0 9 1 0 NaN NaN 1]) newline ...
         'DMIG*   ' pad(seID,16,'right') sprintf('%-16u',[1 0]) newline ];
    A = strrep(A,'NaN','   ');

% Add DMIG vector and remove final \n
    A = [A sprintf(['*' repmat(' ',1,7) '%-16u%-16u %.1E\n'], tmp')];

% Specify double precision
%    A = strrep(A,'E+','D+');

