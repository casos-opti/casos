% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function x = vector_to_indeterminates(obj)
% Convert vector of monomials to indeterminate variables.

[tf,I] = is_monom(obj);

assert(tf, 'Pattern must be vector of monomials.')
assert(obj.mindeg == 1 && obj.maxdeg == 1, 'Pattern must be of degree one.')

% get indeterminate variables
vars = indeterminates(obj);
% return variables in order
x = vars(I);

end
