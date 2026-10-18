% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = nabla(a,x)
% Evaluate polynomial nabla operator (Jacobian matrix).

assert(~is_operator(a), 'Not allowed for operators.')
assert(is_indet(x), 'Second argument must be vector of indeterminates.')

if (isempty(a) || isempty(x))
    % empty Jacobian matrix
    b = a.empty(numel(a),length(x));
    return
end

% else:
b = a.new_poly;

% compute coefficient matrix of Jacobian
[S,b.coeffs] = coeff_nabla(a.get_sparsity,a.coeffs,x);

% set sparsity
b = set_sparsity(b,S);

end
