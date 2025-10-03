% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = sparsify(obj)
% Sparsify polynomial expressions.

b = obj.new_poly;

% sparsify coefficients
coeffs = sparsify(obj.coeffs);

% remove zero terms
[S,b.coeffs] = coeff_update(obj.get_sparsity,coeffs);

% set sparsity
b = set_sparsity(b,S);

end
