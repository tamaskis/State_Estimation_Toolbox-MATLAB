%==========================================================================
%
% obsv_LTI  Observability analysis assuming a linear time-invariant system.
%
%   rank_Ob = obsv_LTI(Fk,Hk)
%
% Author: Tamas Kis
% Last Update: 2022-03-31
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   Fk      - (n×n×N double) discrete dynamics Jacobian at each sample
%   Hk      - (p×p×N double) discrete measurement Jacobian at each sample
%
% -------
% OUTPUT:
% -------
%   rank_Ob - (N×1 double) rank of the observability matrix at each sample
%
% -----
% NOTE:
% -----
%   --> N = number of samples
%   --> Since no measurement is made at the first sample time, we set the
%       rank of the observability matrix at the first sample equal to NaN.
%   --> Since we neglect the control input is not known at the final sample
%       time (we only need the control input at the second to last sample
%       time in order to form the estimate at the last sample time), we set
%       the rank of the observability matrix at the last sample equal to
%       NaN.
%
%==========================================================================
function rank_Ob = obsv_LTI(Fk,Hk)
    
    % number of samples
    N = size(Fk,3);
    
    % preallocates array
    rank_Ob = NaN(N,1);

    % loops over all sample times
    for k = 2:(N-1)
        rank_Ob(k) = rank(obsv(Fk(:,:,k),Hk(:,:,k)));
    end
    
end