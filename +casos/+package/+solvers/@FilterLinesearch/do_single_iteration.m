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
[x_star,dual_star,sol_iter,sol_qp,feas_res_flag,stats] = solve_q_sdp(obj,...
                                                                    iter,...
                                                                    x_k,...
                                                                    p0,...
                                                                    Bk,...
                                                                    args,...
                                                                    stats);

% initilizing timing of backstepping components in case we have to switch
% to feasibility restoration
stats{iter}.time_stats.hess_approx_time    = 0;
stats{iter}.time_stats.total_backstep_time = 0;
stats{iter}.time_stats.filter_accept_time  = 0;
stats{iter}.time_stats.suff_decrease_time  = 0;
stats{iter}.time_stats.soc_time            = 0;


if feas_res_flag
    return % leave single iteration and invoke feasibility restoration
end

%% backtracking filter-linesearch
alpha_max   = obj.opts.filter_struct.alpha_max;

gamma_theta = obj.opts.filter_struct.gamma_theta ;
gamma_phi   = obj.opts.filter_struct.gamma_phi ;

lagrange_filter = obj.opts.filter_struct.lagrange_filter;

alpha = alpha_max;

f_type = 0;
armijo  = 0;


t_backstep = tic;
measCorrectionTime = [];
measSuffDecreaseTime = 0;

counter_linesearch_calls = 0;
counter_correction_calls = 0;

% just for first iteration
alpha_min = 0;

% invoke feasibility restoration if step-length is below minimum
feas_res_flag = 2;

% Backtracking linesearch
while alpha > alpha_min 
  
    counter_linesearch_calls = counter_linesearch_calls + 1;
    % compute search direction (primal and dual)
    dk   = (x_star    - x_k);
    dkl  = (dual_star-dual_k);


    t_filter_accept = tic;
    % check filter acceptance
    [x_k1,theta_x_k1,f_x_k1 ,filter_acceptance,constraint_violation] = check_filter_acceptance(obj,...
                                                                                                linesearch_filter,...
                                                                                                alpha,...
                                                                                                x_k,...
                                                                                                dual_k,...
                                                                                                dk,...
                                                                                                dkl,...
                                                                                                p0,...
                                                                                                args);

    measFilterAcceptTime = toc(t_filter_accept);

    dual_k1 = dual_k + alpha*(dual_star-dual_k);

    % heuristic: there might be blocking entries in the filter:
    % check first-order optimality conditions to potentially aboard
    % linesearch early as possible if we are "optimal"
    if full(obj.eval_gradLag(x_k1,p0,dual_k1))  <= obj.opts.tolerance_opt*max(1,full(obj.eval_gradLag(x_k1,p0,dual_k1))) ...
                                                    && theta_x_k1 <= obj.opts.tolerance_con
        feas_res_flag = 0;
        break
    end

    if lagrange_filter
        L_k1= full(obj.L(x_k1,p0,dual_k1));
        L_k = full(obj.L(x_k,p0,dual_k));
    else
        L_k1 = [];
        L_k  = [];
    end

    if filter_acceptance

        t_suff_decrease = tic;
        % check sufficient decrease
        [suff_decrease_flag,f_type,armijo,linesearch_filter] = check_suff_decrease(obj,...
                                                                                    alpha,...
                                                                                    x_star,...
                                                                                    x_k,...
                                                                                    p0,...
                                                                                    theta_xk,...
                                                                                    theta_x_k1, ...
                                                                                    f_x_k1,...
                                                                                    f_xk,...
                                                                                    L_k1,...
                                                                                    L_k,...
                                                                                    dual_k,...
                                                                                    linesearch_filter);

        measSuffDecreaseTime = toc(t_suff_decrease);

        % either sufficient decrease in cost or constraint violation
        if suff_decrease_flag
            % accept and leave while loop
            feas_res_flag = 0;
            break
        end

        % invoke soc to avoid maratos effect, if necessary
        if ~suff_decrease_flag && alpha == 1  && theta_x_k1 >= theta_xk && obj.opts.enable_soc

            t_soc = tic;
            % compute corrected search direction, compute new trial point and check for filter acceptance
            [x_k1,suff_decrease_soc_flag,f_type,armijo] = second_order_correction(obj,...
                                                                                    x_k,...
                                                                                    x_star,...
                                                                                    p0,...
                                                                                    Bk,...
                                                                                    args,...
                                                                                    linesearch_filter,...
                                                                                    alpha,...
                                                                                    theta_xk,...
                                                                                    f_xk,...
                                                                                    dkl,...
                                                                                    L_k1,...
                                                                                    L_k,...
                                                                                    dual_k);

            measCorrectionTime = toc(t_soc);
            % leave while loop if corrected step is acceptable to filter
            if suff_decrease_soc_flag
                feas_res_flag = 0;
                break % means we found a solution with SOC; leave loop
            end

        end

    else % not acceptable to linesearch_filter

        % invoke soc to avoid maratos effect, if necessary
        if ~filter_acceptance && alpha == 1 && theta_x_k1 >= theta_xk && obj.opts.enable_soc
            counter_correction_calls = counter_correction_calls + 1;
            t_soc = tic;
            % compute corrected search direction, compute new trial point and check for filter acceptance
            [x_k1,suff_decrease_soc_flag,f_type,armijo] = second_order_correction(obj,...
                                                                                    x_k,...
                                                                                    x_star,...
                                                                                    p0,...
                                                                                    Bk,...
                                                                                    args,...
                                                                                    linesearch_filter,...
                                                                                    alpha,...
                                                                                    theta_xk,...
                                                                                    f_xk,...
                                                                                    dkl,...
                                                                                    L_k1,...
                                                                                    L_k,...
                                                                                    dual_k);
                                                                                    
            measCorrectionTime = toc(t_soc);
            if suff_decrease_soc_flag
                feas_res_flag = 0;
                break % means we found a solution with SOC; leave loop
            end

        end

    end

    % if not acceptable to filter and/or no sufficient progress update alpha
    alpha = 1/2*alpha;

    % compute alpha_min
    alpha_min = compute_alpha_min(obj,x_k,dk,p0,theta_xk,linesearch_filter);

    % if alpha < alpha_min  % just avoid very small step lengths
    %     % invoke feasibility restoration if step-length is below minimum
    %     feas_res_flag = 2;
    %     break
    % end

    if alpha_min == 0 && ~filter_acceptance
        % invoke feasibility restoration if step-length is below minimum
        feas_res_flag = 2;
        break
    end

end % end of backtracking linesearch
measTotalBackstepTime = toc(t_backstep);


% only augment if either f-type or armijo are not fulfilled
if (~f_type || ~armijo)

    if lagrange_filter
        % augment filter with a small margin
        linesearch_filter = [linesearch_filter;[theta_xk*(1-gamma_theta), L_k - gamma_phi*theta_xk]];
    else
        % augment filter with a small margin
        linesearch_filter = [linesearch_filter;[theta_xk*(1-gamma_theta), f_xk - gamma_phi*theta_xk]];
    end

end

% output of current iteration
sol_iter.x_k1       = x_k1;
sol_iter.dual_k1    = dual_k1;
sol_iter.theta_x_k1 = theta_x_k1;
sol_iter.f_x_k1     = f_x_k1;
sol_iter.alpha_k    = alpha;
sol_iter.dual_qp    = dual_star;

% store constraint violations for each constraint
stats{iter}.constraint_violation = constraint_violation;

t_hess_approx = tic;
% update quasi-Newton/exact Hessian
if strcmpi(obj.opts.hessian_approx,'BFGS')
    % damped BFGS
    Bk = damped_bfgs(obj,Bk,x_k,p0,sol_iter,iter);

elseif strcmpi(obj.opts.hessian_approx,'levenberg')
    % see Betts book
    H  = full(obj.eval_Hessian(x_k1,p0,dual_k1));
    Bk = casos.package.solvers.SequentialCommon.hessian_levenberg(H);

elseif strcmpi(obj.opts.hessian_approx,'regularization')
    % own regularization method
    H = full(obj.eval_Hessian(x_k1,p0,dual_k1));
    Bk = casos.package.solvers.SequentialCommon.regularize_hessian(H);

elseif strcmpi(obj.opts.hessian_approx,'mirroring')
    % Verschueren
    H = full(obj.eval_Hessian(x_k1,p0,dual_k1));
    [L,D] = ldl(H);

    D(D<0) = 1e-6;
    Bk = L*abs(D)*L';

else
    error('Unknown option "%s" for Hessian approximation.',obj.opts.hessian_approx)
end

measHessApproxTime = toc(t_hess_approx);

% get timing stats
stats{iter}.time_stats.hess_approx_time    = measHessApproxTime;
stats{iter}.time_stats.total_backstep_time = measTotalBackstepTime;
stats{iter}.time_stats.filter_accept_time  = measFilterAcceptTime;
stats{iter}.time_stats.suff_decrease_time  = measSuffDecreaseTime;
stats{iter}.counter_correction_calls       = counter_correction_calls;
stats{iter}.counter_linesearch_calls       = counter_linesearch_calls;

if ~isempty(measCorrectionTime)
    stats{iter}.time_stats.soc_time        = measCorrectionTime;
else
    stats{iter}.time_stats.soc_time        = 0;
end


end
