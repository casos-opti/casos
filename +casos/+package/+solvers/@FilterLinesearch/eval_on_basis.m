% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function argout = eval_on_basis(obj,argin)
% Solve sequential SOS problem.
t_solve = tic;
import casos.package.UnifiedReturnStatus

%% print output with current problem and setting
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug',['\t\t Ca' char(931) 'oS-Nonlinear sum-of-squares optimization suite \n']);
printf(obj.log,'debug','\t\t Sequential Quadratic Sum-of-Squares Solver v1.1, October''26 \n');
printf(obj.log,'debug','\t\t GNU GENERAL PUBLIC LICENSE \n');
printf(obj.log,'debug','\t\t Institute of Flight Mechanics and Control, Univ. of Stuttgart \n');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','\n');
printf(obj.log,'debug','Problem\n ')
printf(obj.log,'debug','Decision variables \t\t= %d (coefficents) \n',obj.display_para.no_decVar);
printf(obj.log,'debug',' SOS constraints \t\t= %d \n', obj.display_para.no_sosCon);
printf(obj.log,'debug',' Linear constraints \t\t= %d\n',0 );

printf(obj.log,'debug','-------------------------------------------------------\n')
printf(obj.log,'debug','Settings \n ')
printf(obj.log,'debug','Solver\t\t\t\t- %s\n ', obj.display_para.solver);
printf(obj.log,'debug','Maximum iterations\t\t- %d\n ', obj.opts.max_iter);
printf(obj.log,'debug','Constraint violation Check \t- %s\n ',obj.opts.con_vio_check);
printf(obj.log,'debug','Second-order correction \t- %s \n',mat2str(obj.opts.enable_soc));
printf(obj.log,'debug',' Hessian-approximation\t\t- %s\n',obj.opts.hessian_approx);
printf(obj.log,'debug',[' ' char(949) '_con\t\t\t\t- %d\n'],obj.opts.tolerance_con);
printf(obj.log,'debug',[' ' char(949) '_opt\t\t\t\t- %d\n'],obj.opts.tolerance_opt);
printf(obj.log,'debug','-------------------------------------------------------\n');
printf(obj.log,'debug','\n');
printf(obj.log,'debug','Solver build time \t\t:  %0.2d s\n', obj.display_para.solver_build_time);
printf(obj.log,'debug','Start solving ... \n');
printf(obj.log,'debug','\n');

printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');



t_init = tic;
% input arguments
args = argin;

% initial guess (user)
x_k = args{1};

% parameter from nonlinear problem
p0   = args{2};

% initialize iteration statistics struct; + 1 because we add also the
% zeroth iteration/initial step
stats.single_iterations = cell(1,obj.opts.max_iter+1);

stats.options = obj.opts;

% initialze dual variables
dual_k = zeros(obj.init_para.no_dual_var,1);

% select hessian approximation method
% Note: If BFGS is used we do not pre-compute the Hessian to save build
% time! Hence, we use the identy matrix and potentially user scaling
if strcmpi(obj.opts.hessian_approx,'BFGS')
        Bk =  eye(obj.init_para.size_B)*obj.opts.scale_bfgs0;
else
    
    % if other than BFGS are used we either use identity with scaling or
    % pre-computed hessian and regularize it
    if strcmpi(obj.opts.hessian_init,'identity')
        % initialize BFGS-matrix = identity x scaling
        Bk =  eye(obj.init_para.size_B)*obj.opts.scale_bfgs0;
    
    elseif strcmpi(obj.opts.hessian_init,'analytical')
        % initialize Hessian approximation with regularization
        H = full(obj.eval_Hessian(x_k,p0,dual_k));
    
        Bk = casos.package.solvers.SequentialCommon.regularize_hessian(H);
    
    else
        error('Unknown option "%s" for Hessian initialization.',obj.opts.hessian_init)
    end
end


if strcmpi(obj.opts.con_vio_check,'signed-distance')
    % initialize filter and first iterate
    args_con_vio    =  args;
    args_con_vio{2} =  [p0; x_k];
    args_con_vio{3} = -inf(obj.init_para.con_vio.no_con,1);
    args_con_vio{4} =  inf(obj.init_para.con_vio.no_con,1);

    args_con_vio{6} = [];
    args_con_vio{7} = [];

    % compute constraint violation of initial guess
    sol_con_vio = eval_on_basis(obj.solver_con_vio, args_con_vio);

    lb = args{6};
    ub = args{7};
    lin_violation = obj.linvio(x_k,p0,lb,ub);

    % largest constraint violation
    theta_x_k = full(max(0,max([sol_con_vio{1}; lin_violation])));

elseif strcmpi(obj.opts.con_vio_check,'projection')
    % initialize filter and first iterate
    args_con_vio    =  args;
    args_con_vio{2} =  [p0; x_k];
    args_con_vio{3} = [];
    args_con_vio{4} =  [];

    args_con_vio{6} = [];
    args_con_vio{7} = [];

    % compute constraint violation of initial guess
    sol_con_vio = eval_on_basis(obj.solver_con_vio, args_con_vio);

    lb = args{6};
    ub = args{7};
    lin_violation = obj.linvio(x_k,p0,lb,ub);

    % largest constraint violation: cost gives us minimum. violation
    theta_x_k = full(max(0,max([sol_con_vio{2}; lin_violation])));

else
    g_val = full(obj.eval_constraint_samples(x_k,obj.opts.user_sample));

    min_vio = min(min(g_val));
    if min_vio < 0
        theta_x_k = abs(min_vio);
    else
        theta_x_k = 0;
    end
end

% evaluate cost function of nonlinear problem
f_x_k     = full(obj.eval_cost(x_k,p0));

% add some boundary to it
linesearch_filter = [max(1,theta_x_k)*10, inf];

iter = 1;

% just for display output in first iteration
[delta_xi_double,~,delta_dual_double,~] = obj.eval_step_norms(x_k,x_k,dual_k,dual_k);
delta_xi_double   = full(delta_xi_double);
delta_dual_double = full(delta_dual_double);
alpha_k           = 1;

% counter until solution is marked almost optimal
counter_almost_opt_iter    = 0;

% flag to indicate current iterate is feasible ( 1 = feasible)
feasibility_flag          = 1;

% counter how often feasibility restoration is called
counter_restoration_calls = 0;



% first entry is the initialization step i.e. at initialization 
% no line search
stats.single_iterations{1}.counter_correction_calls        = 0;
stats.single_iterations{1}.counter_linesearch_calls        = 0;

% in case of signed-distance we provide the statistics of the underlying
% SDP and to find out which constraint(s) might be the blocking point
if strcmpi(obj.opts.con_vio_check,'signed-distance')

    stats.single_iterations{1}.sdp_sol                         = obj.solver_con_vio.get_stats;
    stats.single_iterations{1}.constraint_violation            = [sol_con_vio{1}; lin_violation];

else
    stats.single_iterations{1}.sdp_sol                         = [];
    stats.single_iterations{1}.constraint_violation            = [];
end

stats.single_iterations{1}.time_stats.main_sdp_time        = 0;
stats.single_iterations{1}.time_stats.hess_approx_time     = 0;
stats.single_iterations{1}.time_stats.total_backstep_time  = 0;
stats.single_iterations{1}.time_stats.filter_accept_time   = 0;
stats.single_iterations{1}.time_stats.suff_decrease_time   = 0;
stats.single_iterations{1}.time_stats.soc_time             = 0;
stats.single_iterations{1}.time_stats.feas_res_time        = 0;


% end of time measurement for initialization
stats.init_time = toc(t_init);


% main-loop
while iter <= obj.opts.max_iter

    % display output header every 10th iteration
    if ~mod(iter,10) && iter > 0
        printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    end

    % evaluate first-order optimality and termination thresholds once per iterate
    gradLag_k = full(obj.eval_gradLag(x_k,p0,dual_k));
    [opt_thres,opt_thres_10x,opt_thres_self] = obj.eval_termination(f_x_k, full(obj.eval_gradLag_inf_norm(x_k,p0,dual_k)), gradLag_k, delta_xi_double, obj.opts.tolerance_opt);
    opt_thres      = full(opt_thres);
    opt_thres_10x  = full(opt_thres_10x);
    opt_thres_self = full(opt_thres_self);

    % display output current iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , gradLag_k   );
    
    % store current iterate: cost, max. constraint violation, current step
    % length and first-order optimality condition
    stats.single_iterations{iter}.iter         = iter-1;
    stats.single_iterations{iter}.cost         = f_x_k;
    stats.single_iterations{iter}.max_con_vio  = theta_x_k;
    stats.single_iterations{iter}.step_length  = alpha_k;
    stats.single_iterations{iter}.d_ldx        = gradLag_k;
   
    % store counter of current iterate
    stats.single_iterations{iter}.counter_restoration_calls = counter_restoration_calls; 

    % initialize work statistics of the NEXT entry
    stats.single_iterations{iter+1}.counter_correction_calls        = 0;
    stats.single_iterations{iter+1}.counter_linesearch_calls        = 0;
    stats.single_iterations{iter+1}.sdp_sol                         = [];
    stats.single_iterations{iter+1}.constraint_violation            = [];
    stats.single_iterations{iter+1}.time_stats.main_sdp_time        = 0;
    stats.single_iterations{iter+1}.time_stats.hess_approx_time     = 0;
    stats.single_iterations{iter+1}.time_stats.total_backstep_time  = 0;
    stats.single_iterations{iter+1}.time_stats.filter_accept_time   = 0;
    stats.single_iterations{iter+1}.time_stats.suff_decrease_time   = 0;
    stats.single_iterations{iter+1}.time_stats.soc_time             = 0;
    stats.single_iterations{iter+1}.time_stats.feas_res_time        = 0;


    %% check convergence (first-order optimality)
    if gradLag_k <= opt_thres  ... % scaled optimality
            && theta_x_k <= obj.opts.tolerance_con && feasibility_flag == 0                                                                                     % constraint violation

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found\n');
        measSolveTime = toc(t_solve);
        printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);
        feasibility_flag = 1;
        stats.solution_status = 'Optimal solution';
        obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        break

    end

    % check if current iterate is solved to an acceptable level
    if gradLag_k <= opt_thres_10x  ... % scaled optimality ...
            && theta_x_k <= obj.opts.tolerance_con

        % increase counter
        counter_almost_opt_iter = counter_almost_opt_iter +1;

        % we have 5-10 in a row, seems we solved it to an acceptable level/almost optimal
        if counter_almost_opt_iter >= obj.opts.almost_opt_count
            printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
            printf(obj.log,'debug','Solution status: Almost optimal.\n');
            measSolveTime = toc(t_solve);
            printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
            printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
            printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);

            feasibility_flag = 1;

            stats.solution_status = 'Alomst Optimal solution.';
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
            break
        end
    else
        % if current iterate is not solved to an acceptable level, reset counter
        counter_almost_opt_iter = 0;
    end

    feasibility_flag = 0;

    % perform linesearch 
    % Note: Work statistics are written to entry iter+1, i.e. the entry of the iterate this step produces
    [sol_iter,sol_qp,feas_res_flag,stats.single_iterations,obj,linesearch_filter,Bk] = do_single_iteration(obj, ...
                                                                                                iter+1,...
                                                                                                x_k,...
                                                                                                dual_k,...
                                                                                                theta_x_k,...
                                                                                                f_x_k,...
                                                                                                Bk,...
                                                                                                p0,...
                                                                                                args, ...
                                                                                                linesearch_filter, ...
                                                                                                stats.single_iterations);

    
    stats.single_iterations{iter+1}.time_stats.feas_res_time = 0;

    % Invoke feasibility restoration if necessary
    if feas_res_flag >= 1 && ...                   % restoration requested
            theta_x_k >  obj.opts.tolerance_con    % is the constraint violation larger then tolerance; otherwise restoration does not make sense
        
        counter_restoration_calls = counter_restoration_calls + 1;

        % parameter (just for readabilty)
        gamma_theta = obj.opts.filter_struct.gamma_theta ;
        gamma_phi   = obj.opts.filter_struct.gamma_phi ;

        % augment filter with iterate where restoration is invoked; this point
        % shall be avoided in the future
        linesearch_filter = [linesearch_filter;[theta_x_k*(1-gamma_theta), f_x_k - gamma_phi*theta_x_k]];

        time_stats_backup = stats.single_iterations{iter+1}.time_stats;
        t_feas_res = tic;
        % feasibility restoration phase
        [sol_iter,sol_qp_feas,linesearch_filter,feasibility_flag,stats.single_iterations] = obj.feas_res_solver.eval_extended(linesearch_filter,...
                                                                                                                                x_k,...
                                                                                                                                theta_x_k,...
                                                                                                                                p0,...
                                                                                                                                obj,...
                                                                                                                                stats.single_iterations,...
                                                                                                                                iter+1);

        stats.single_iterations{iter+1}.time_stats = time_stats_backup;
        stats.single_iterations{iter+1}.time_stats.feas_res_time = toc(t_feas_res);


        if feasibility_flag >=1

            % continue from restored solution
            [~,delta_xi_double,delta_dual_double,~] = obj.eval_step_norms(x_k,sol_iter.x_k1,dual_k,dual_k);
            delta_xi_double   = full(delta_xi_double);
            delta_dual_double = full(delta_dual_double);

            x_k                 = sol_iter.x_k1;
            sol_iter.dual_k1    = dual_k;            % we keep dual soultion from previous iterate
            theta_x_k           = sol_iter.theta_x_k1;
            f_x_k               = sol_iter.f_x_k1;
            alpha_k             = sol_iter.alpha_k;

        else
            % could not recover return iterate where restoration was
            % invoked; return point where invoked
            sol    = sol_qp;

            sol{1} = x_k;
            sol{5} = dual_k;

            break
        end

    elseif feas_res_flag >= 1 &&  theta_x_k < obj.opts.tolerance_con
        % we are already better then treshhold i.e. restoration can't
        % make it better
        sol    = sol_qp;

        % if one, means Q-SDP is infeasible, thus sol.iter should be
        % zero hence return "old" iterate
        if feas_res_flag == 1
            sol{1} = x_k;
            sol{5} = dual_k;
        else
            sol{1} = sol_iter.x_k1;
            sol{5} = sol_iter.dual_k1;
        end

        % just a final check
        if gradLag_k <= opt_thres_self ...
                && theta_x_k <= obj.opts.tolerance_con

            printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
            printf(obj.log,'debug','Solution status: Optimal solution found\n');
            measSolveTime = toc(t_solve);
            printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
            printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
            printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);
            feasibility_flag = 1;

        else
            % stalled
            feasibility_flag = -2;

        end

        break



    else % continue with normal iteration

        % for display output
        [~,delta_xi_double,~,delta_dual_double] = obj.eval_step_norms(x_k,sol_iter.x_k1,dual_k,sol_iter.dual_k1);
        delta_xi_double   = full(delta_xi_double);
        delta_dual_double = full(delta_dual_double);

        % set values for next iterate
        x_k       = sol_iter.x_k1;
        dual_k    = sol_iter.dual_k1;
        theta_x_k = sol_iter.theta_x_k1;
        f_x_k     = sol_iter.f_x_k1;
        alpha_k   = sol_iter.alpha_k;

    end

    % store counter in statistics structure
    stats.single_iterations{iter}.counter_restoration_calls = counter_restoration_calls; 
    
    % return last solution
    if ~isempty(sol_qp)
        sol    = sol_qp;
    elseif ~isempty(sol_qp_feas)
        sol = sol_qp_feas;
    else
        sol = [];
    end

    % exchange the primal/dual solution of QP with the accepted iterate
    sol{1} = sol_iter.x_k1;
    sol{2} = sol_iter.f_x_k1;
    % sol{3} = obj.eval_constraint(sol_iter.x_k1,p0);
    % dual variables estimated by underlying quadratic SDP
    sol{4} = sol_qp{4};
    sol{5} = sol_iter.dual_k1;

    % increase iteration counter
    iter = iter + 1;

end % end of while loop

% evaluate first-order optimality and termination threshold once for the final iterate
gradLag_k = full(obj.eval_gradLag(x_k,p0,dual_k));
[~,~,opt_thres_self] = obj.eval_termination(f_x_k, full(obj.eval_gradLag_inf_norm(x_k,p0,dual_k)), gradLag_k, delta_xi_double, obj.opts.tolerance_opt);
opt_thres_self = full(opt_thres_self);

% maximum number of iterations exhausted: store the last iterate (index
% max_iter+1 corresponds to iteration max_iter); its work statistics are
% already stored by the last loop pass
if iter > obj.opts.max_iter
    stats.single_iterations{iter}.iter         = iter-1;
    stats.single_iterations{iter}.cost         = f_x_k;
    stats.single_iterations{iter}.max_con_vio  = theta_x_k;
    stats.single_iterations{iter}.step_length  = alpha_k;
    stats.single_iterations{iter}.d_ldx        = gradLag_k;
    stats.single_iterations{iter}.counter_restoration_calls = counter_restoration_calls;
end

% final display output for user
if iter >= obj.opts.max_iter

    % display output header
    if ~mod(iter,10) && iter > 0
        printf(obj.log,'debug','%-8s%-15s%-15s%-15s%-15s%-10s%-10s\n', 'iter', 'obj', ['||' char(916) 'x'  '||_inf'], ['||' char(916) char(955)  '||_inf'], [char(920) '(x_k)'],char(945),'||dLdx||');
        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    end

    % display output current iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , gradLag_k   );


    if gradLag_k <= opt_thres_self ...
            && theta_x_k <= obj.opts.tolerance_con && feasibility_flag == 0

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Optimal solution found\n');
        measSolveTime = toc(t_solve);
        printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);

        stats.solution_status = 'Optimal Solution';
        obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
       

    else

        printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
        printf(obj.log,'debug','Solution status: Maximum number of iterations reached\n');
        measSolveTime = toc(t_solve);
        printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
        printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
        printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);

        stats.solution_status = 'Max. Iterations';
        obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
       
    end
end

if feasibility_flag == 0 && iter < obj.opts.max_iter

    % print last iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , gradLag_k   );

    printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    printf(obj.log,'debug','Solution status: Problem infeasible.\n');
    measSolveTime = toc(t_solve);
    printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
    printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
    printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);

    stats.solution_status = 'Problem infeasible';
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;

elseif feasibility_flag == -1 && iter < obj.opts.max_iter


    % print last iterate
    printf(obj.log,'debug','%-8d%-15e%-15e%-15e%-15e%-10f%-10e\n',...
        iter-1, f_x_k , delta_xi_double, delta_dual_double, theta_x_k , alpha_k , gradLag_k   );


    printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    printf(obj.log,'debug','Solution status: Solver stalled in feasibility restoration.\n');
    measSolveTime = toc(t_solve);
    printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
    printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
    printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);

    stats.solution_status = 'Solver stalled';
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
    
elseif feasibility_flag == -2 && iter < obj.opts.max_iter

    printf(obj.log,'debug','------------------------------------------------------------------------------------------\n');
    printf(obj.log,'debug','Solution status: Feasible solution. Restoration not invoked because already feasible.\n');
    measSolveTime = toc(t_solve);
    printf(obj.log,'debug',['Solution time: ' num2str(measSolveTime) ' s\n']);
    printf(obj.log,'debug',['Build time: ' num2str(obj.display_para.solver_build_time) ' s\n']);
    printf(obj.log,'debug',['Total time: ' num2str(obj.display_para.solver_build_time+measSolveTime) ' s\n']);


    stats.solution_status = 'Feasible solution';
    obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
end

stats.single_iterations = stats.single_iterations(1:iter);

stats.total_solve_time    = measSolveTime;
stats.solver_build_time   = obj.display_para.solver_build_time;
stats.iterations        = iter-1; % we do not count the initial guess
stats.time_per_iteration  = stats.total_solve_time/stats.iterations;

% assign output
argout = sol;

% store iteration statistics
obj.sequential_stats = stats;

end % end of function
