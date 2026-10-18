% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_plus(obj,S2,coeff1,coeff2)
% Add two operator coefficient matrices.

% expand internally
[cfa,cfb,Si,So] = expand_internal(obj,S2,coeff1,coeff2);

% remove zero terms
[coeffs,Si,So] = removeZero(cfa+cfb,Si,So);

% new sparsity pattern
S = new_from_coefficients(coeffs,Si,So);

end
