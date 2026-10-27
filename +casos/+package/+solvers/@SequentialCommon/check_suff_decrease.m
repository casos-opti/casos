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
s_theta     = obj.opts.filter_struct.s_theta ;
s_phi       = obj.opts.filter_struct.s_phi ;
gamma_theta = obj.opts.filter_struct.gamma_theta ;
gamma_phi   = obj.opts.filter_struct.gamma_phi ;


delta           = obj.opts.filter_struct.delta;
eta             = obj.opts.filter_struct.eta;
lagrange_filter = obj.opts.filter_struct.lagrange_filter;

theta_min = min(linesearch_filter(1,1),1)*obj.opts.filter_struct.theta_min;

% search direction
dk = x_star    - x_k;

% descent direction
nabla_f_dir = full(obj.eval_nabla_f_dir(x_k,x_star,p0));

% f-type switching: descent direction + progress in cost better than con. violation
f_type = alpha*nabla_f_dir  < 0 && (-alpha*nabla_f_dir)^s_phi*alpha^(1-s_phi) > delta*theta_xk^s_theta;

% initialize armijo boolean
armijo  = 0;

if  f_type && theta_xk <= theta_min

    % check armijo-condition
    if lagrange_filter
        armijo = L_k1 <= L_k +eta*full(( obj.eval_dLdx(x_k,p0,dual_k))'*dk);
    else
        armijo = f_x_k1 <= f_xk + eta*alpha*nabla_f_dir;
    end

    if armijo
        suff_decrease_flag = 1;
    else
        suff_decrease_flag = 0;
    end

else

    % check progress w.r.t. previous iteration
    if theta_x_k1 <= (1-gamma_theta)*theta_xk  || f_x_k1 <= f_xk - gamma_phi*theta_xk      %  L_k1 <= L_k - gamma_phi*theta_xk
        suff_decrease_flag = 1;
    else
        suff_decrease_flag = 0;
    end

end % end of progress check

end % end of function