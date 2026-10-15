% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_update(obj,coeffs,sz,dim)
% Update polynomial coefficient matrix.

if nargin > 3
    % finish after matrix operation
    coeffs = finishMatrixOp(coeffs,sz,dim);
elseif nargin < 3
    % matrix dimensions do not change
    sz = size(obj);
end

assert(all(size(coeffs) == [obj.nterm prod(sz)]),'Notify the developers.')

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,obj.degmat,obj.indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,sz);

end
