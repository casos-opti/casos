% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = evaluate(a,b)
% Evaluate operator.

assert(is_operator(a), 'Evaluation only allowed for operators.')

assert(all(size(a.sparsity_in) == size(b)), 'Dimension mismatch for operator evaluation.')

c = a.new_poly;

% compute dot product
[S,c.coeffs] = coeff_evaluate(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
