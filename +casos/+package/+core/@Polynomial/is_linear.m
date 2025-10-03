% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function tf = is_linear(p,q)
% Check if polynomial is linear in symbols.

assert(is_symbolic(q),'Second argument must be purely symbolic.')

% get nonzero coordinates
Q = coordinates(q);

% check if coefficients are linear
tf = is_linear(p.coeffs,Q);

end
