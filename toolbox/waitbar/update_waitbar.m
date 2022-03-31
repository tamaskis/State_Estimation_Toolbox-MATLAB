%==========================================================================
%
% update_waitbar  Updates the waitbar.
%
%   prop = update_waitbar(i,N,wb,prop)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   wb      - (1×1 Figure) waitbar
%   n       - (1×1 double) current sample number (i.e. iteration)
%  	N       - (1×1 double) total number of samples (i.e. iterations)
%   prop    - (1×1 double) cutoff proportion to trigger waitbar update
%
% -------
% OUTPUT:
% -------
%   prop    - (1×1 double) cutoff proportion to trigger waitbar update
%
% -----
% NOTE:
% -----
%   --> "prop" is an integer multiple of 0.1 so that the waitbar is
%       only updated after every additional 10% of progress.
%
%==========================================================================
function prop = update_waitbar(i,N,wb,prop)
    
    % only updates waitbar if current proportion exceeds cutoff prop.
    if i/N > prop
        
        % updates waitbar
        waitbar(i/N,wb);
        
        % updates cutoff proportion needed to trigger waitbar update
        prop = prop+0.1;
        
    end
    
end