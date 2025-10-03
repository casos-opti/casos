% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function p = new_sym(obj,dstr,varargin)
% Create a new symbolic polynomial (helper function).

if nargin < 3
    % constant scalar
    S = casos.Sparsity.scalar;

elseif ~isa(varargin{1},'casos.Sparsity')
    % zero-degree symbol
    S = casos.Sparsity(varargin{:});

elseif nargin > 3
    % repeat scalar pattern
    assert(isscalar(varargin{1}),'Sparsity pattern must be scalar.')

    S = casos.Sparsity.dense(varargin{:});

else
    % sparsity pattern
    S = varargin{1};
end

p = obj.new_poly;

% assign sparsity
p = set_sparsity(p,S);
% symbolic coefficient matrix
p.coeffs = sym_coeff(obj,dstr,coeff_sparsity(S));

end
