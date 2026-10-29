% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularizes a given Hessian matrix to ensure it is positive definite.
% The implemented method is inspired by the Levenberg method described in
% Section 2.5 in "Betts, J.- Practical Methods for Optimal Control and
% Estimation Using Nonlinear Programming, 2010". A small modification is
% done by interpolating tau parameter.
%
% Inputs:
%   H   - (n x n) Hessian matrix (must be symmetric)
%   tau - Regularization parameter (small positive scalar, e.g., 1e-6)
%
% Output:
%   H_reg - Regularized positive definite Hessian matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H_regularized = hessian_levenberg(H)

epsilon = sqrt(eps);

% Diagonal elements of H
h_diag = diag(H);

% Sum of absolute values of elements in each row
row_sum = sum(abs(H), 2);

% Gershgorin bound (vectorized)
gershgorin_bound = h_diag - (row_sum - abs(h_diag));

% Take the minimum Gershgorin bound
sigma = min(gershgorin_bound);

% Determine tau based on sigma
tau = max(0, epsilon - sigma) / (abs(sigma) + 1);


% Add the regularization term to the Hessian
I = eye(size(H,1)); % Identity matrix
H_regularized = H + tau * (abs(sigma) + 1) * I;
end
