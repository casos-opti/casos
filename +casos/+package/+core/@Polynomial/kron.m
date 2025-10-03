% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = kron(a,b)
% Compute Kronecker product of two polynomial matrices.

assert(~is_operator(a) && ~is_operator(b), 'Not allowed for operators.')

% handle simple case(s) for speed up
if isempty(a) || isempty(b)
    % product with empty polynomial is empty

    c = a.zeros(size(a).*size(b));
    return
end

% else
c = a.new_poly;

% Kronecker product of coefficients
[S,c.coeffs] = coeff_kron(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

% set sparsity
c = set_sparsity(c,S);

end
