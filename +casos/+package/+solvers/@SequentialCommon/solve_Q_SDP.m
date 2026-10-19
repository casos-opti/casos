% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description:
%
% Call the actual quadratic subproblem (SDP with quadratic cost) of the
% form
%   min 1/2*d'*B*d + ∇f'*d
%   st. g(x) + ∇g(x)⋅d in Σ[x]
%       x in R[x]
%
% where d = (x^* - x_k) is the search direction, B is the approximation of
% the Hessian, ∇(⋅) the first-derivative and (⋅)^* is the optimal solution.
%
% Remark: All decision variables are setup as regular polynomials. SOS
% polynomials are enforced via constraints!
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_star,dual_star,sol_iter,sol_qp,feas_res_flag,stats] = solve_Q_SDP(obj, ...
    iter, ...
    x_k, ...
    p0, ...
    Bk, ...
    args, ...
    stats)

import casos.package.UnifiedReturnStatus

% pre-pare input for Q-SDP
args{2}  = [p0; x_k; Bk(:)];
% adjust bounds
% args{3}  = args{3} - x_k;
% args{4}  = args{4} - x_k;
measSolver =tic;
sol_qp   = eval_on_basis(obj.solver_convex, args);
solveTimeSDP = toc(measSolver);
% store iteration info
stats{iter}.sdpsol           = obj.solver_convex.get_stats  ;
stats{iter}.timeStats.mainSDP = solveTimeSDP;

%% check solution status
switch (obj.solver_convex.get_stats.UNIFIED_RETURN_STATUS)

    % optimal solution found (or at least almost feasible)
    case UnifiedReturnStatus.SOLVER_RET_SUCCESS

        x_star    = sol_qp{1};
        dual_star = sol_qp{5};
        feas_res_flag = 0;
        sol_iter = [];

        % case where the solver could not decide or numerical issues
    case UnifiedReturnStatus.SOLVER_RET_LIMITED

        feas_res_flag = 1;
        x_star        = x_k;
        dual_star     = [];
        sol_iter      = [];

        % strictly infeasible
    otherwise
        feas_res_flag = 1;
        x_star    = x_k;
        dual_star = [];
        sol_iter = [];

end

end