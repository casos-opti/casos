% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function S = tuplet(n,m,i,j,x,degrees)
% Create sparsity pattern with tuplets.

assert(n >= 0 && m >= 0,'Dimensions must be nonnegative.')
assert(length(i) == length(j),'Mismatch of length of indices.')
assert(all(i >= 0) && all(j >= 0),'Indices must be nonnegative.')
assert(all(i < n) && all(j < m),'Indices out of range (%dx%d).',n,m)
assert(is_indet(x),'x must be indeterminate variables.')
assert(length(x) == size(degrees,2),'Mismatch of number of indeterminate variables.')
assert(length(i) == size(degrees,1),'Mismatch of length of degrees.')

S = casos.Sparsity;

% convert subindices to column index
col = sub2ind([n m],i+1,j+1); % 1-index

% coefficients with nonunique degree matrix
coeffs = casadi.Sparsity.triplet(length(col),n*m,0:length(col)-1,col-1); % 0-index

% make degree matrix unique
[coeffs,degmat] = uniqueDeg(coeffs,degrees);

% remove redundant zeros
[S.coeffs,S.degmat,S.indets] = removeZero(coeffs,degmat,casos.Indeterminates(x));

% set dimensions
S.matdim = [n m];

end
