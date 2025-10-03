% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function coeffs = sparsity_prod(coeffs,sz,dim)
% Compute product of coefficient sparsity pattern along given dimension.

cfa = pattern_inverse(coeffs);

% compute sum of inverse
cfb = sparsity_sum(cfa,sz,dim);

% product is inverse of sum
coeffs = pattern_inverse(cfb);

end
