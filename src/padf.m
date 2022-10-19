% Faster than pad, can also truncate the char array // same in xbdfread
function str = padf(str,w)
    if ~exist('w','var')    w = 72;                 end                     % Set desired char array width
    if ~ischar(str)         str = char(str);        end                     % Convert to char array

	if size(str,2)==w                                                       % Do nothing
        return
    elseif size(str,2)<w                                                    % Expand to 72 fields per line, if less
        if      ndims(str)==2 	str(:,size(str,2)+1:w) = ' ';               % Faster than (:,~,:)
        elseif  ndims(str)==3 	str(:,size(str,2)+1:w,:) = ' ';
        end
    else                                                                    % Remove anything after char 72, if lines are longer
        if      ndims(str)==2	str(:,w+1:end) = '';
        elseif  ndims(str)==3   str(:,w+1:end,:) = '';
        end
	end
end

