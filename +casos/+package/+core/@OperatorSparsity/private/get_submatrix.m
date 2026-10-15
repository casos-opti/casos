% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function M = get_submatrix(coeffs,I,J)
% Return the IxJ submatrix of the coefficients matrix.

if isa(coeffs,'casadi.Sparsity')
    % subreference of sparsity pattern
    M = sub(coeffs,I-1,J-1);    % CasADi uses 0-index

else
    % subreference coefficient matrix
    M = coeffs(I,J);
end

end
