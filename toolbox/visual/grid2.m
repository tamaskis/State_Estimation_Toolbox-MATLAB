%==========================================================================
%
% grid2  TODO.
%
%   TODO
%
% Author: Tamas Kis
% Last Update: 2021-12-09
%
% REFERENCES:
%   [1] https://www.mathworks.com/matlabcentral/fileexchange/36953-patchline
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   t           - (1×N double) time vector
%   x           - (1×N double) a posteriori state variable estimates
%   x_lower     - (1×N double) lower bound on state variable estimates
%   x_upper     - (1×N double) upper bound on state variable estimates
%   pp          - (struct) plot parameters
%
%==========================================================================
function grid2
    
    % axis handle 
    H = gca;

    % axis limits
    xmin = H.XLim(1);
    xmax = H.XLim(2);
    ymin = H.YLim(1);
    ymax = H.YLim(2);

    % vectors of tick marks
    xtick = H.XTick;
    ytick = H.YTick;

    % coordinates for vertical grid lines
    vx = zeros(1,3*length(xtick));
    vy = zeros(size(vx));
    for i = 1:length(xtick)
        j = 3*(i-1)+1;
        vx(j) = xtick(i);
        vy(j) = ymin;
        vx(j+1) = xtick(i);
        vy(j+1) = ymax;
        vx(j+2) = NaN;
        vy(j+2) = NaN;
    end

    % coordinates for horizontal grid lines
    hx = zeros(1,3*length(xtick));
    hy = zeros(size(hx));
    for i = 1:length(ytick)
        j = 3*(i-1)+1;
        hx(j) = xmin;
        hy(j) = ytick(i);
        hx(j+1) = xmax;
        hy(j+1) = ytick(i);
        hx(j+2) = NaN;
        hy(j+2) = NaN;
    end

    % plots grid lines
    hold on;
    patch([hx(:);NaN],[hy(:);NaN],[0.15,0.15,0.15],'edgealpha',0.15);
    patch([vx(:);NaN],[vy(:);NaN],[0.15,0.15,0.15],'edgealpha',0.15);

end