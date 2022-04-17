%==========================================================================
%
% transform_affine  Affine transformation of Gaussian random vectors.
%
%   [y,Py] = transform_affine(x,P,A)
%   [y,Py] = transform_affine(x,P,A,b)
%
% Author: Tamas Kis
% Last Update: 2022-04-01
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (n×N double) a posteriori state estimates
%   P       - (n×n×N double) a posteriori error covariances
%   A       - (1×1 function_handle) linear portion of affine 
%             transformation, Aₖ = A(xₖ,k) (A : ℝⁿ×ℤ → ℝᵐˣⁿ)
%   b       - (1×1 function_handle) (OPTIONAL) translation portion of
%             affine transformation, bₖ = A(xₖ,k) (b : ℝⁿ×ℤ → ℝᵐ) (defaults
%             to 0)
%
% -------
% OUTPUT:
% -------
%   y       - (m×N double) a posteriori alternate state estimates
%   Py      - (m×m×N double) a posteriori alternate error covariances
%
%==========================================================================
function [y,Py] = transform_affine(x,P,A,b)
    
    % alternate state dimension
    m = size(A(x(:,1),0),1);    

    % defaults "b" to 0 if not input
    if (nargin < 4) || isempty(b)
        b = @(xk,k) zeros(m,1);
    end

    % number of samples
    N = size(x,2);

    % preallocates arrays
    y = zeros(m,N);
    Py = zeros(m,m,N);

    % alternate state estimates at every sample
    for k = 0:(N-1)

        % index for accessing arrays (switch from 0- to 1-based indexing)
        kk = k+1;

        % parameters defining affine transformation at kth sample time
        Ak = A(x(:,kk),k);
        bk = b(x(:,kk),k);

        % alternate state estimates at kth sample time
        y(:,kk) = Ak*x(:,kk)+bk;
        Py(:,:,kk) = Ak*P(:,:,kk)*Ak.';

    end

end