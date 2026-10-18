% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,I] = op_subsref(obj,ii)
% Subreference into nonzero coefficients.

S = coeff_subsref(obj,obj.coeffs,ii,[length(ii) 1]);

% get indices of referenced nonzeros
[~,jc] = get_triplet(obj.coeffs);
% find referenced coefficients
I = find(ismember(jc+1,ii));    % CasADi has 0-index

end
