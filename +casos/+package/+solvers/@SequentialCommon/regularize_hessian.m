% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularizes a given Hessian matrix to ensure it is positive definite.
% The implemented method regularizes the Hessian by computung the frobenius
% norm and to ensure positive definitness. A "new" Hessian is than
% constructed.
%
% Inputs:
%   H   - (n x n) Hessian matrix (must be symmetric)
%   tau - Regularization parameter (small positive scalar, e.g., 1e-6)
%
% Output:
%   H_reg - Regularized positive definite Hessian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_reg = regularize_hessian(H)

    if nargin < 2
        epsilon = sqrt(eps);          % relative floor, ~1.5e-8
    end

    H = (H + H') / 2;                 % enforce symmetry

    [Q, lambda] = eig(H, 'vector');   % orthonormal Q for symmetric H
    tau = epsilon * max(norm(H, 'fro'), 1);

    if min(lambda) < tau
        lambda_reg = max(lambda, tau);

        H_reg = Q * diag(lambda_reg) * Q';   % nearest PD in Frobenius norm

        H_reg = (H_reg + H_reg') / 2;

    else

        H_reg = H;
        
    end

end
