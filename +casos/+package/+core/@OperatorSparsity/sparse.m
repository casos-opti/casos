% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/casos-opti/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function S = sparse(varargin)
% Create all-sparse operator pattern.

assert(nargin > 0, 'Not enough input arguments.')
assert(nargin < 3, 'Too many input arguments.')

if isnumeric(varargin{1})
    % input-output dimensions given
    Si = casos.Sparsity(varargin{end},1);
    So = casos.Sparsity(varargin{1},1);

else
    % input-output patterns given
    Si = casos.Sparsity(varargin{1});
    So = casos.Sparsity(varargin{end});
end

% coefficients are sparse
M = casadi.Sparsity(nnz(So),nnz(Si));

% create operator pattern
S = casos.package.core.OperatorSparsity(M,Si,So);

end
