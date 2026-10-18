% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function S = to_vector(obj,I,row_vector)
% Convert a scalar sparsity pattern into vector of terms.

assert(isscalar(obj),'Dimensions must be 1x1.')

nt = obj.nterm;

if nargin < 2 || isempty(I)
    % default ordering
    I = 1:nt;
end
if nargin < 3 || ~row_vector
    % return column vector
    sz = [nt 1];
else
    % return row vector
    sz = [1 nt];
end

% new sparsity pattern
coeffs = casadi.Sparsity.triplet(nt,nt,I-1,(1:nt)-1); % NOTE: Casadi has zero-based index
S = new_from_coefficients(coeffs,obj.degmat,obj.indets,sz);

end
