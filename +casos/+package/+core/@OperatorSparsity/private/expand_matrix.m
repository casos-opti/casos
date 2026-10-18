% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function M = expand_matrix(S_coeffs,coeffs,sz,Iin,Iout)
% Expand operator matrix to joint nonzeros.

[ii,jj] = get_triplet(S_coeffs);

% new sparsity pattern of operator matrix
SpA = casadi.Sparsity.triplet(sz(1),sz(2),Iout(ii+1)-1,Iin(jj+1)-1);

% return expanded operator matrix
M = sparsity_cast(coeffs,SpA);

end
