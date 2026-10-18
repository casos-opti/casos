% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,I] = op_subsref(obj,ii)
% Subreference into nonzero coefficients.

[S,~,I1,I2] = coeff_subsref(obj,obj.sparsity_M,ii,[]);

% get indices of referenced nonzeros
[ic,jc] = get_triplet(obj.sparsity_M);
% find referenced coefficients
I = find(ismember(ic+1,I1) & ismember(jc+1,I2));    % CasADi has 0-index

end
