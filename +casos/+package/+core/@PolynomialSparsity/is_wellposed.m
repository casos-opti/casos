% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function tf = is_wellposed(obj)
% Check if sparsity is well posed.

tf = size(obj.coeffs,2) == prod(obj.matdim)      ... number of elements
    && size(obj.coeffs,1) == size(obj.degmat,1)  ... number of terms
    && size(obj.degmat,2) == length(obj.indets); ... number of variables

end
