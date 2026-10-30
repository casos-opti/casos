% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description:
%
% Check if new iterate is acceptable to filter. Compute constraint
% violation and cost at new iterate. Compare the pair {cost constraint
% vio.} to filter entries. If at least one of them is better in cost or
% constrain violation, than accept.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [x_k1,theta_x_k1,f_x_k1 ,filter_acceptance,all_violations] = check_filter_acceptance(obj, ...
    linesearch_filter, ...
    alpha, ...
    x_k, ...
    dual_k, ...
    dk, ...
    dkl, ...
    p0, ...
    args)

% use Langrange function instead of cost function in filter
lagrange_filter = obj.opts.filter_struct.lagrange_filter;

if strcmpi(obj.opts.con_vio_check,'signed-distance')

    % violation of the linear constraints
    lb = args{6}; 
    ub = args{7}; 

    % compute new solution candidate
    if lagrange_filter
        [x_k1,lin_violation,L_k1] = obj.eval_x_k1(x_k,p0,alpha,dk,dkl,dual_k,lb,ub);
        x_k1 = full(x_k1);
        lin_violation = full(lin_violation);
        L_k1 = full(L_k1);
    else
        [x_k1,lin_violation,f_x_k1] = obj.eval_x_k1(x_k,p0,alpha,dk,lb,ub);
        x_k1 = full(x_k1);
        f_x_k1 = full(f_x_k1);
    end

    % compute constraint violation of new solution candidate
    args_con_vio    =  args;
    args_con_vio{2} =  [p0; x_k1];
    args_con_vio{3} = -inf(obj.init_para.con_vio.no_con,1);
    args_con_vio{4} =  inf(obj.init_para.con_vio.no_con,1);

    args_con_vio{6} = [];
    args_con_vio{7} = [];

    sol_con_vio = eval_on_basis(obj.solver_con_vio, args_con_vio);

    % extract signed-distances
    all_violations = [sol_con_vio{1}; lin_violation];

    % get the largest one
    theta_x_k1 = full(max(0,max(all_violations)));
    
elseif strcmpi(obj.opts.con_vio_check,'projection')

    % violation of the linear constraints
    lb = args{6}; 
    ub = args{7}; 

    % compute new solution candidate
    if lagrange_filter
        [x_k1,lin_violation,L_k1] = obj.eval_x_k1(x_k,p0,alpha,dk,dkl,dual_k,lb,ub);
        x_k1 = full(x_k1);
        lin_violation = full(lin_violation);
        L_k1 = full(L_k1);
    else
        [x_k1,lin_violation,f_x_k1] = obj.eval_x_k1(x_k,p0,alpha,dk,lb,ub);
        x_k1 = full(x_k1);
        f_x_k1 = full(f_x_k1);
    end

    % compute constraint violation of new solution candidate
    args_con_vio    =  args;
    args_con_vio{2} =  [p0; x_k1];
    args_con_vio{3} =  []; % projection has no constraint
    args_con_vio{4} =  [];

    args_con_vio{6} = [];
    args_con_vio{7} = [];

    sol_con_vio = eval_on_basis(obj.solver_con_vio, args_con_vio);

    % extract projection: cost gives us minimum. violation
    all_violations = [sol_con_vio{2}; lin_violation];

    % get the largest one
    theta_x_k1 = full(max(0,max(all_violations)));

else
    % evalaute at current solution and samples
    if lagrange_filter
        [x_k1,all_violations,theta_x_k1,L_k1] = obj.eval_x_k1(x_k,p0,alpha,dk,dkl,dual_k,obj.opts.user_sample);
        x_k1 = full(x_k1);
        all_violations = full(all_violations);
        theta_x_k1 = full(theta_x_k1);
        L_k1 = full(L_k1);
    else
        [x_k1,all_violations,theta_x_k1,f_x_k1] = obj.eval_x_k1(x_k,p0,alpha,dk,obj.opts.user_sample);
        x_k1 = full(x_k1);
        all_violations = full(all_violations);
        theta_x_k1 = full(theta_x_k1);
        f_x_k1 = full(f_x_k1);
    end
end

% check filter acceptance
% get filter inputs
theta_l = linesearch_filter(:,1);
f_l     = linesearch_filter(:,2);

% new point lies in forbidden region if both are larger than filter entries
dominance_bool      = [];

% boundary belongs to forbidden region i.e. >=
dominance_bool(:,1) = theta_x_k1 >= theta_l;

if lagrange_filter
    dominance_bool(:,2) = L_k1     >= f_l;
else
    dominance_bool(:,2) = f_x_k1  >= f_l;
end

% check pairs; if both one, means pair lies in forbidden region
dominance_bool = all(dominance_bool, 2);

if any(dominance_bool) % means not acceptable to linesearch_filter
    filter_acceptance = 0;
else
    filter_acceptance = 1;
end

end
