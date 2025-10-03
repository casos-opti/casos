% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function out = str_terms(obj)
% Print element-wise monomial terms.

[~,L] = get_degmat(obj);
monom = str_monoms(obj);

% prepare output
out = cell(size(obj));
% preassign chars
out(:) = {'()'};

for i = find(obj)
    % build element-wise monomial terms
    m = join(monom(L(i,:)),',');

    if sum(L(i,:)) > 1
        % embed vector in parentheses
        out(i) = compose('(%s)',m{:});
    else
        out(i) = m;
    end
end

end
