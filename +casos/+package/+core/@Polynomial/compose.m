% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = compose(a,b)
% Composition of operators.

assert(is_operator(a) && is_operator(b), 'Composition only allowed for operators.')

assert(all(size(a.sparsity_in) == size(b.sparsity_out)), 'Dimension mismatch for operator composition.')

c = a.new_poly;

% compute dot product
[S,c.coeffs] = coeff_compose(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
