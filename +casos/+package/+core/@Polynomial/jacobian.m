% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function b = jacobian(a,x)
% Compute symbolic Jacobian matrix of vector polynomial expression.
%
% Use of JACOBIAN for the polynomial Jacobian matrix is not allowed. 
% Use NABLA instead.

if is_indet(x)
    % not allowed: compute polynomial Jacobian matrix.
    error('Use NABLA for the polynomial Jacobian.')
end

% else:
assert(is_symbolic(x),'Second argument must be symbolic polynomial.')

% project x to basis
[X,Si] = coordinates(x);

% project f to basis
[A,So] = coordinates(a);

% G is jacobian of coefficients
B = jacobian(A,X);

% return operator
b = a.new_poly(B,Si,So);

end
