% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs,I1,I2] = coeff_subsref(obj,coeffs,ii,~)
% Subreference operator coefficient matrix.

[ia,ja] = ind2sub(size(obj),ii);

% convert subindices into linear indices of input/output patterns
Ia = unique(ia);
Ja = unique(ja);

assert(numel(ii) == length(Ia)*length(Ja),'Size of subreference must not change.')

% subreference input/output patterns
[Si,I2] = op_subsref(obj.sparsity_in,Ja);
[So,I1] = op_subsref(obj.sparsity_out,Ia);

if isa(coeffs,'casadi.Sparsity')
    % subreference sparsity pattern
    coeffs = sub(coeffs,I1,I2,true);

else
    % subreference coefficient matrix
    coeffs = coeffs(I1,I2);
end

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
