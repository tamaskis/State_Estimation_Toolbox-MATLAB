%==========================================================================
%
% initialize_waitbar  Initializes the waitbar.
%
%   prop = initialize_waitbar(i,N,wb,prop)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   wb_in   - (char or 1×1 logical) (OPTIONAL) waitbar parameters
%               --> input as "true" if you want waitbar with default 
%                   message displayed
%               --> input as a char array storing a message if you want a
%                   custom message displayed on the waitbar
%   msg     - (char) default messagae (if not specified by wb_in)
%   m       - (1×1 double) parameter number for function that is calling
%             "initialize_waitbar" that controls 
%   n       - (1×1 double) current sample number (i.e. iteration)
%  	N       - (1×1 double) total number of samples (i.e. iterations)
%   prop    - (1×1 double) cutoff proportion to trigger waitbar update
%
% -------
% OUTPUT:
% -------
%   wb      - (1×1 Figure) waitbar
%   prop    - (1×1 double) cutoff proportion to trigger waitbar update
%   disp_wb - (1×1 logical) true if waitbar is turned on, false otherwise
%
% -----
% NOTE:
% -----
%   --> If "disp_wb" is false, then "wb" and "prop" are returned as NaN.
%
%==========================================================================
function [wb,prop,display_waitbar] = initialize_waitbar(wb_in,msg)
    
    % determines if waitbar is on or off
    if (islogical(wb_in) && ~wb_in)
        display_waitbar = false;
    else
        display_waitbar = true;
    end

    % resets the waitbar message if specified by "wb_in"
    if display_waitbar && ischar(wb_in)
        msg = wb_in;
    end

    % initialize cutoff proportion needed to trigger waitbar update to 0.1
    prop = NaN;
    if display_waitbar
        prop = 0.1;
    end

    % initializes the waitbar
    wb = NaN;
    if display_waitbar
        wb = waitbar(0,msg);
    end
    
end