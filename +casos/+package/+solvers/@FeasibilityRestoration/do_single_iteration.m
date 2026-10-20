% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [sol_iter,sol_qp,feas_res_flag,stats,obj,linesearch_filter,Bk] = do_single_iteration(obj, ...
    iter,...
    x_k,...
    dual_k,...
    theta_xk,...
    f_xk,...
    Bk,...
    p0,...
    args, ...
    linesearch_filter, ...
    stats)

%% evaluate convex SOS problem and check feasibility
[x_star,dual_star,sol_iter,sol_qp,feas_res_flag,stats] = solve_q_sdp(obj,iter,x_k,p0,Bk,args,stats);


if feas_res_flag
    return % leave single iteration and invoke feasibility restoration
end

%% backtracking filter-linesearch
alpha_max   = obj.opts.filter_struct.alpha_max;

gamma_theta = obj.opts.filter_struct.gamma_theta ;
gamma_phi   = obj.opts.filter_struct.gamma_phi ;

lagrange_filter = obj.opts.filter_struct.lagrange_filter;


alpha = alpha_max;

% Backtracking linesearch
while true

    % compute search direction
    dk   = (x_star    - x_k);
    dkl  = (dual_star-dual_k);

    % check filter acceptance
    [x_k1,theta_x_k1,f_x_k1 ,filter_acceptance,constraint_violation] = check_filter_acceptance(obj,linesearch_filter,alpha,x_k,dual_k,dk,dkl,p0,args);

    dual_k1 = dual_k + alpha*(dual_star-dual_k);

    if lagrange_filter
        L_k1= full(obj.L(x_k1,p0,dual_k1));
        L_k = full(obj.L(x_k,p0,dual_k));
    else
        L_k1 = [];
        L_k  = [];
    end

    if filter_acceptance

        % check sufficient decrease
        [suff_decrease_flag,f_type,armijo,linesearch_filter] = check_suff_decrease(obj,alpha,x_star,x_k,p0,theta_xk,theta_x_k1, f_x_k1,f_xk,L_k1,L_k,dual_k,linesearch_filter);

        if suff_decrease_flag
            % accept and leave while loop
            break
        end

        % invoke soc to avoid maratos effect, if necessary
        if ~suff_decrease_flag && alpha == 1  && theta_x_k1 >= theta_xk

            % compute corrected search direction, compute new trial point and check for filter acceptance
            [x_k1,suff_decrease_flag,f_type,armijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,linesearch_filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k);

            % leave while loop if corrected step is acceptable to filter
            if suff_decrease_flag
                break % means we found a solution with SOC; leave loop
            end

        end


    else % not acceptable to linesearch_filter

        % invoke soc to avoid maratos effect, if necessary
        if ~filter_acceptance && alpha == 1 && theta_x_k1 >= theta_xk

            % compute corrected search direction, compute new trial point and check for filter acceptance
            [x_k1,suff_decrease_flag,f_type,armijo] = second_order_correction(obj,x_k,x_star,p0,Bk,args,linesearch_filter,alpha,theta_xk,f_xk,dkl,L_k1,L_k,dual_k);

            if suff_decrease_flag
                break % means we found a solution with SOC; leave loop
            end

        end

    end

    % if not acceptable to filter and/or no sufficient progress update alpha
    alpha = 1/2*alpha;

    % compute alpha_min
    alpha_min = compute_alpha_min(obj,x_k,dk,p0,theta_xk,linesearch_filter);

    if alpha < alpha_min  || alpha <= eps
        % invoke feasibility restoration if step-length is below minimum
        feas_res_flag = 2;
        break
    end


end % end of backtracking linesearch


% only augment if new iterate is accepted and either f-type or armijo are not fulfilled
if suff_decrease_flag && (~f_type || ~armijo)

    if lagrange_filter
        % augment filter with a small margin
        linesearch_filter = [linesearch_filter;[theta_xk*(1-gamma_theta), L_k - gamma_phi*theta_xk]];
    else
        % augment filter with a small margin
        linesearch_filter = [linesearch_filter;[theta_xk*(1-gamma_theta), f_xk - gamma_phi*theta_xk]];
    end

end


% output of current iterate
sol_iter.x_k1       = x_k1;
sol_iter.dual_k1    = dual_k1;
sol_iter.theta_x_k1 = theta_x_k1;
sol_iter.f_x_k1     = f_x_k1;
sol_iter.alpha_k    = alpha;
sol_iter.dual_qp    = dual_star;
stats{iter}.constraint_violation = constraint_violation;

% we only support BFGS because we have acutally a linear cost + a quadratic
% term; other methods faild in the past
Bk = damped_bfgs(obj,Bk,x_k,p0,sol_iter,iter);



end
