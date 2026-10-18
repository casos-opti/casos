% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_subsref(obj,coeffs,ii,sz)
% Subreference polynomial coefficient matrix.

nt = obj.nterm;

if isa(coeffs,'casadi.Sparsity')
    % subreference sparsity pattern
    coeffs = sub(coeffs,1:nt,ii,true);

else
    % subreference coefficient matrix
    coeffs = coeffs(:,ii);
end

% remove coefficients, degrees, and/or indeterminates 
% that do not appear in the referenced polynomial
[coeffs,degmat,indets] = removeZero(coeffs,obj.degmat,obj.indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,sz);

end
