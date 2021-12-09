%==========================================================================
%
% process_noise_covariance  Obtain the process noise covariance.
%
%   Q = process_noise_covariance(x_true,fd)
%
% Author: Tamas Kis
% Last Update: 2021-11-15
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x_true  - (n×T) ground truth state vector time history
%   f       - (function_handle) f(x,u,k) --> (f:Rn×Rm×R->Rn) discrete-time 
%                                            nonlinear dynamics equation
%   u       - (m×T double) control input time history
%
% -------
% OUTPUT:
% -------
%   Q       - (n×n double) process noise covariance
%
%==========================================================================
function Q = process_noise_covariance(x_true,f,u)
    
    dx_true = diff(x_true')';

    % preallocates array for simulated differences
    dx_sim = zeros(size(dx_true));

    % TODO
    for k = 1:size(dx_sim,2)
        dx_sim(:,k) = f(x_true(:,k),u(:,k),k);
    end

    difference = dx_true-dx_sim;

    Q_vec = zeros(size(x_true,1),1);

    for i = 1:length(Q_vec)
        Q_vec(i) = std(difference(i,:));
    end
    
    Q_vec
    
    Q = diag(Q_vec);
    

end