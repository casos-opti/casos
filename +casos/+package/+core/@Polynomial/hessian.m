% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = hessian(a,x)
% Compute symbolic Hessian matrix of vector polynomial expression.


if is_indet(x)
    % not allowed: compute polynomial Hessian matrix.
    error('Use NABLA for the polynomial Jacobian.')
end

assert(is_symbolic(x),'Second argument must be symbolic polynomial.')
assert(isscalar(a) && is_zerodegree(a),'Objective must be scalar variable.')

% project x to basis
[X,zi] = poly2basis(x);

% project f to basis
A = poly2basis(a);

% G is jacobian of coefficients
B = hessian(A,X);

% return operator
b = a.new_poly(B,zi,zi);

end
