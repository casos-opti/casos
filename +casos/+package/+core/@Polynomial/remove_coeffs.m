% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = remove_coeffs(obj,tol)
% Remove coefficients that are below a given tolerance.

b = obj.new_poly;

% identify coefficients to keep
S_coeffs = sparsity(sparsify(abs(obj.coeffs) >= tol));

% project onto coefficient sparsity
coeffs = project(obj.coeffs,S_coeffs);

% update coefficients
[S,b.coeffs] = coeff_update(obj.get_sparsity,coeffs);

% set sparsity
b = set_sparsity(b,S);

end
