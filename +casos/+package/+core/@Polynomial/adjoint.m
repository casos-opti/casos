% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = adjoint(a)
% Return adjoint operator.

b = a.new_poly;

% adjoint coefficient matrix
[S,b.coeffs] = coeff_adjoint(a.get_sparsity,a.coeffs);

% set sparsity
b = set_sparsity(b,S);

end
