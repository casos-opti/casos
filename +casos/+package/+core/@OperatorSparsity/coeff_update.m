% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_update(obj,coeffs,sz,dim) %#ok<INUSD>
% Update operator coefficient matrix.

if nargin > 2
    % finish after matrix operation
    error('Not supported.')
end

assert(all(size(coeffs) == [nnz(obj.sparsity_out) nnz(obj.sparsity_in)]),'Notify the developers.')

% TODO: remove zero terms
[coeffs,Si,So] = removeZero(coeffs,obj.sparsity_in,obj.sparsity_out);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
