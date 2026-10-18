% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function coeffs = evalMatrixOp(coeffs,op,dim)
% Evaluate matrix operation along specified dimension.

switch (dim)
    case 'all'
        % matrix operation over all elements
        coeffs = op(coeffs,2);

    case {1 2}
        % matrix operation along dimension
        coeffs = op(coeffs,dim);

    otherwise
        error('Invalid dimension input.')
end

end
