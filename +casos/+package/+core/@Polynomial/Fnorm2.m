% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function r = Fnorm2(obj)
% Compute the squared Frobenius norm of a polynomial operator.

assert(is_operator(obj), 'Only allowed for operators.')

r = dot(obj.coeffs,obj.coeffs);

end
