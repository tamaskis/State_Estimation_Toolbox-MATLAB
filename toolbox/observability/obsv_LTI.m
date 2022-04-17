%==========================================================================
%
% obsv_LTI  Observability analysis assuming a linear time-invariant system.
%
%   [rank_Ob,obsvbl] = obsv_LTI(F,H)
%
% Author: Tamas Kis
% Last Update: 2022-04-02
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   F       - (n×n double) discrete dynamics Jacobian
%   H       - (p×n double) discrete measurement Jacobian
%
% -------
% OUTPUT:
% -------
%   rank_Ob - (1×1 double) rank of the observability matrix
%   obsvbl  - (1×1 logical) true if system is observable, false otherwise
%
%==========================================================================
function [rank_Ob,obsvbl] = obsv_LTI(F,H)

    % rank of the observability matrix
    rank_Ob = rank(obsv(F,H));

    % state dimension
    n = size(F,1);

    % determines if system is observable
    if rank_Ob == n
        obsvbl = true;
    else
        obsvbl = false;
    end

end