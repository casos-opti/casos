% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_setnz(obj,I,v)
% Set nonzero elements of a operator coefficients matrix.

sz = coeff_size(obj);
[ii,jj] = coeff_triplet(obj);

% new coefficient sparsity
S_coeffs = casadi.Sparsity.triplet(sz(1),sz(2),ii(I),jj(I));

if nargin > 2
    % cast onto nonzeros
    coeffs = sparsity_cast(v,S_coeffs);
else
    % nothing to do
    coeffs = S_coeffs;
end

% remove zero terms
[coeffs,Si,So] = removeZero(coeffs,obj.sparsity_in,obj.sparsity_out);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
