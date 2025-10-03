% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function C = times_coeffs(A,B)
% Multiply coefficient matrices.

if isa(A,'casadi.Sparsity')
    % intersect sparsity
    C = intersect(A,B);

else
    % multiply coefficients element-wise
    C = times(A,B);
end

end
