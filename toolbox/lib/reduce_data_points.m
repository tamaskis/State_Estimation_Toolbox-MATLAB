%==========================================================================
%
% reduce_data_points  Reduces the number of data points in a data set to a
% specified number.
%
%   [x1_new,...,xn_new] = reduce_data_points(x1,...,xn,N)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-04-16
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   [x1,...,xn]         - original data set with N₀ data points
%                           --> each vector is a 1×N₀ or N₀×1 double
%   N                   - (1×1 double) desired number of data points
%
% -------
% OUTPUT:
% -------
%   [x1_new,...,xn_new] - updated data set with N data points
%                           --> each vector is a 1×N or N×1 double
%
% -----
% NOTE:
% -----
%   --> n = dimension of a data point
%           • For example, if a single data point is represented by an
%             ordered triple (x1,x2,x3), then n = 3.
%   --> N₀ = original number of data points
%   --> N = (desired) new number of data points
%   --> Sometimes, the function will not be able to return exactly N points
%       (due to rounding issues).
%   --> The main purpose of this function is to reduce the size of a data 
%       set when not all the points are needed. For example, plotting 
%       y = x^2 with 100 points rather than with 1000 points will (to the 
%       naked eye) be visually identical, but will be a lot faster for the 
%       computer to perform.
%
%==========================================================================
function varargout = reduce_data_points(varargin)
    
    % extracts data set and number of data points to save
    data = varargin(1:(end-1));
    N = varargin{end};
    
    % determines original number of data points
    N0 = length(data{1});
    
    % determines indices of data points to save
    save_ratio = round(N0/N);
    save = 1:save_ratio:N0;
    
    % preallocates output argument
    varargout = cell(1,nargin-1);
    
    % shrinks data set
    for i = 1:(nargin-1)
        varargout{i} = data{i}(save);
    end
    
end