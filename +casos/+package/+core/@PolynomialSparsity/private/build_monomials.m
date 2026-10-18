% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function z = build_monomials(degmat,indets)
% Build a monomial sparsity pattern.

% number of monomials
nt = size(degmat,1);

% set output
z = casos.package.core.PolynomialSparsity( ...
                           casadi.Sparsity.dense(nt,1), ...
                           sparse(degmat), ...
                           indets, ...
                           [1 1] ...
);

end
