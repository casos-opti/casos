% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_adjoint(obj,coeffs)
% Return coefficient matrix of adjoint.

S = dualize(obj);

% return coefficient matrix for dualization
coeffs = sparsity_cast(coeffs,coeff_sparsity(S));

end
