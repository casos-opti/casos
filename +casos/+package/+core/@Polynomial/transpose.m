% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = transpose(a)
% Transpose of polynomial matrix.

assert(~is_operator(a), 'Not allowed for operators. Use "adjoint" instead.')

b = a.new_poly;

% transpose coefficient matrix
[S,b.coeffs] = coeff_transpose(a.get_sparsity,a.coeffs);

% set sparsity
b = set_sparsity(b,S);

end
