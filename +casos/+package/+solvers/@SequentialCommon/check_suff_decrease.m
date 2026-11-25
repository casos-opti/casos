% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description:
%
% Check if new iterate fulfills sufficient decrease condition in cost and
% ensure that progress (decrease) in constraint violation is smaller than
% cost (f-type switching) and if Amijo condition is
% fulfilled. If not f-type, check if at least sufficient progress
% (both cost or constraint violation) with respect to previous iterate is
% given. If neither f-type nor progress with respect to previous
% iterate, than there is no sufficient progress. In that case either a
% second-order correction is invoked or the step-length must be reduced.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [suff_decrease_flag,f_type,armijo,linesearch_filter] = check_suff_decrease(obj, ...
    alpha, ...
    x_star, ...
    x_k, ...
    p0, ...
    theta_xk, ...
    theta_x_k1, ...
    f_x_k1, ...
    f_xk, ...
    L_k1, ...
    L_k, ...
    dual_k, ...
    linesearch_filter)

% depack options just for readibility
lagrange_filter = obj.opts.filter_struct.lagrange_filter;

theta_min = min(linesearch_filter(1,1),1)*obj.opts.filter_struct.theta_min;

% descent direction
nabla_f_dir = full(obj.eval_nabla_f_dir(x_k,x_star,p0));

if lagrange_filter
    [f_type,armijo_cond,suff_decrease_cond] = obj.eval_suff_decrease(alpha,nabla_f_dir,theta_xk,theta_x_k1,f_x_k1,f_xk,L_k1,L_k,x_star,dual_k,x_k,p0);
else
    [f_type,armijo_cond,suff_decrease_cond] = obj.eval_suff_decrease(alpha,nabla_f_dir,theta_xk,theta_x_k1,f_x_k1,f_xk);
end

f_type = full(f_type);
armijo_cond = full(armijo_cond);
suff_decrease_cond = full(suff_decrease_cond);

% initialize armijo boolean
armijo  = 0;

if  f_type && theta_xk <= theta_min

    % check armijo-condition
    armijo = armijo_cond;

    if armijo
        suff_decrease_flag = 1;
    else
        suff_decrease_flag = 0;
    end

else

    % check progress w.r.t. previous iteration
    if suff_decrease_cond
        suff_decrease_flag = 1;
    else
        suff_decrease_flag = 0;
    end

end % end of progress check

end % end of function
