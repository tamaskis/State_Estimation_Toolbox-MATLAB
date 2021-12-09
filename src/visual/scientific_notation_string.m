%==========================================================================
%
% scientific_notation_string  Returns a string of a number written using
% scientific notation.
%
%   str = scientific_notation_string(x)
%   str = scientific_notation_string(x,n)
%
% Author: Tamas Kis
% Last Update: 2021-12-09
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (1×1 double) number
%   n       - (OPTIONAL) (1×1 double) number of digits after decimal point 
%             (i.e. precision)
%
% -------
% OUTPUT:
% -------
%   str     - (string) string storing number in scientific notation, 
%             formatted to be used with LaTeX interpreter (but does NOT 
%             include leading and trailing "$")
%   
%==========================================================================
function str = scientific_notation_string(x,n)

    % defaults "n" to 4
    if (nargin < 2), n = 4; end
    
    % convert number to scientific notation string with n decimal places
    string = sprintf("%."+n+"e",x);
    
    % split the string where "e" is
    string_parts = strsplit(string,'e');

    % isolates the coefficient and exponent
    coefficient = str2double(string_parts(1));
    exponent = str2double(string_parts(2));

    % return the string (formatted for use with LaTeX interpreter)
    str = sprintf("%."+n+"f",coefficient)+"\times10^{"+exponent+"}";
    
end