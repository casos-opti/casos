% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = sum(a,dim)
% Return sum of array elements.

assert(~is_operator(a), 'Not allowed for operators.')

if nargin > 1
    % nothing to do
elseif isempty(a) && ~isvector(a)
    % sum of empty matrix (return scalar)
    dim = 'all';
elseif isrow(a)
    % sum along second dimension (return column)
    dim = 2;
else
    % sum along first dimension (return row)
    dim = 1;
end

if isempty(a)
    % sum of empty matrix is zero
    b = a.zeros(sizeofMatrixOp(a.get_sparsity,dim));
    return
end

% else
b = a.new_poly;

% compute sum of coefficients
[S,b.coeffs] = coeff_sum(a.get_sparsity,a.coeffs,dim);

% set sparsity
b = set_sparsity(b,S);

end
