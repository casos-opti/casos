% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function S = new_from_coefficients(coeffs,Si,So)
% Create sparsity pattern with (generic) coefficients.

if ~isa(coeffs,'casadi.Sparsity')
    coeffs = sparsity(coeffs);
end

% new sparsity pattern
S = casos.Sparsity.create(casos.package.core.OperatorSparsity(coeffs,Si,So));

end
