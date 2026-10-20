% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function argout = eval_on_basis(obj,argin)

args = argin;

xi_k = args{1};

% parameter from nonlinear problem
p0   = args{2};


Bk = eye(obj.size_hessian(1));

% initialize iteration info
stats = cell(1,obj.opts.max_iter);


args{2}  = [p0; xi_k; Bk(:);zeros(length(xi_k),1)];
% evaluate convex SOS problem
sol = eval_on_basis(obj.sossolver, args);

% store iteration info
stats{1} = obj.sossolver.get_stats;


% return last solution
argout = sol;

% store iteration info
obj.sequential_stats.iter = stats;

end
