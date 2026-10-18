% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function v = coeff_getnz(obj,coeffs,I)
% Return nonzero elements of an operator coefficients matrix.

if isa(coeffs,'casadi.Sparsity')
    % return dense sparsity pattern
    v = casadi.Sparsity.dense(numel(I));
    return
end

% else:
nz = sparsity_cast(coeffs,casadi.Sparsity.dense(nnz(obj),1));

% select nonzeros
v = nz(I);

end
