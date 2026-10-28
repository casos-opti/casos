% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (Abstract) SequentialCommon < casos.package.solvers.SosoptCommon
    % Base class for sequential sum-of-squares algorithms.

    properties (Constant,Access=protected)
        SequentialOptions = [casos.package.solvers.SosoptCommon.sosopt_options
            {'sossol', 'The convex sum-of-squares solver to be used in the subproblem.'
            'sossol_options', 'Options to be passed to the SOS solver.'
            'tolerance_con', 'Absolute tolerance for stopping criterion of constraint violation.'
            'tolerance_opt', 'Absolute tolerance for stopping criterion of constraint violation.'
            'filter_struct', 'Structure containing parameter for filter linesearch.'
            'scale_bfgs0','Scaling parameter for initial BFGS matrix.'
            'enable_soc', 'Flag to turn on off the second-order-correction [default: true]'
            'hessian_init','Method to initialize Hessian.'
            'hessian_approx','Hessian (Lagrangian) approximation method'
            'max_iter', 'Maximum number of iterations.'
            'almost_opt_count','Number of iterations to check for almost optimal. '
            'feasibility_restoration','Control cost function of feasibility restoration.'
            'con_vio_check','Select method for constraint violation check'
            'user_sample','If sampling is used, user can provide samples.'
            'verbose', 'Turn on/off iteration display.'}
            ];

        
    end

    properties (Access=protected)
        % low-level solvers
        solver_convex;
        solver_con_vio;
        solver_soc;

        % evaluation of linear cone constraints
        linvio;

        % Lagrangian and derivative
        eval_L;
        eval_dLdx;

        % linesearch
        eval_gradCost;
        eval_constraint_samples;
        eval_nabla_f_dir;

        % damped BFGS
        damped_bfgs;
        eval_s;
        eval_y;
        eval_r;

        % convergence check
        eval_cost;
        eval_constraint;
        eval_gradLag;
        eval_gradLag_inf_norm;
        eval_Hessian;
        eval_step_norms;
        eval_termination;

        % parameters
        init_para;
        feas_res_para;
        display_para;

        % logging & info
        log;
        status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        sequential_stats = struct('iter',[]);
    end

    methods (Static)
        function options = get_options
            % Return static options.
            options = casos.package.solvers.SequentialCommon.SequentialOptions;
        end
    end

    methods (Access=protected)
        % prepare problem structures
        buildproblem(obj,nlsos);
        
        % iteration for overloading
        varargout = run_iteration(varargin);
        varargout = do_single_iteration(varargin);

        % check if the sufficient decrease condition for the filter are met
        varargout = check_suff_decrease(varargin);

        % internal evaluation
        argout = eval_on_basis(obj,argin);
    end

    methods (Static)
        % iteration for overloading/hessian approximation
        varargout = regularize_hessian(varargin);
        varargout = hessian_levenberg(varargin);
        varargout = hessian_regularization(varargin);
    end

    methods
        % Constructor
        function obj = SequentialCommon(name,nlsos,varargin)
            obj@casos.package.solvers.SosoptCommon(name,nlsos,varargin{:});

            % states
            if ~isfield(nlsos,'x')
                error('No decision variables given.')
            else
                nlsos.x = casos.PS(nlsos.x);
            end
            % parameter
            if ~isfield(nlsos,'p')
                % not parametrized
                nlsos.p = casos.PS;
            else
                nlsos.p = casos.PS(nlsos.p);
            end
            % objective
            if ~isfield(nlsos,'f')
                nlsos.f = casos.PS(0);
            else
                nlsos.f = casos.PS(nlsos.f);
            end
            % constraints
            if ~isfield(nlsos,'g')
                nlsos.g = casos.PS;
            else
                nlsos.g = casos.PS(nlsos.g);
            end

            % make sure cholseky option is passed
            % sossol_options0.sdpsol_options.cholesky_method = 'num';

            % ------------------  default options ------------------
            % select underlying conic solver
            if ~isfield(obj.opts,'sossol'), obj.opts.sossol                                     = 'mosek'; end

            % options for the underlying convex problem
            if ~isfield(obj.opts,'sossol_options'), obj.opts.sossol_options                     = struct; end

            % tolerance for first order KKT conditions 
            if ~isfield(obj.opts,'tolerance_con'), obj.opts.tolerance_con                       = 1e-6; end
            if ~isfield(obj.opts,'tolerance_opt'), obj.opts.tolerance_opt                       = 1e-4; end

            % --- Hessian approximation ---
            % select method
            if ~isfield(obj.opts,'hessian_approx'), obj.opts.hessian_approx                     = 'regularization'; end
            % how to initialize hessian
            if ~isfield(obj.opts,'hessian_init'),   obj.opts.hessian_init                       = 'analytical'; end
            %initial scaling factor for Hessian 
            if ~isfield(obj.opts,'scale_bfgs0'), obj.opts.scale_bfgs0                           = 1; end

            % feasibility restoration
            if ~isfield(obj.opts,'feasibility_restoration'), obj.opts.feasibility_restoration   = 'regularize'; end

            % flag to turn on/off second order correction
            if ~isfield(obj.opts,'enable_soc'), obj.opts.enable_soc                             = true; end

            % structure with parameter for filter
            if ~isfield(obj.opts,'filter_struct'), obj.opts.filter_struct                             = struct; end
            if ~isfield(obj.opts.filter_struct,'alpha_max'),       obj.opts.filter_struct.alpha_max       = 1; end
            if ~isfield(obj.opts.filter_struct,'s_theta'),         obj.opts.filter_struct.s_theta         = 1.1; end
            if ~isfield(obj.opts.filter_struct,'s_phi'),           obj.opts.filter_struct.s_phi           = 2.3; end
            if ~isfield(obj.opts.filter_struct,'gamma_theta'),     obj.opts.filter_struct.gamma_theta     = 1e-5; end
            if ~isfield(obj.opts.filter_struct,'gamma_phi'),       obj.opts.filter_struct.gamma_phi        = 1e-5; end
            if ~isfield(obj.opts.filter_struct,'gamma_alpha'),     obj.opts.filter_struct.gamma_alpha     = 0.05; end
            if ~isfield(obj.opts.filter_struct,'delta'),           obj.opts.filter_struct.delta           = 0.8; end
            if ~isfield(obj.opts.filter_struct,'eta'),             obj.opts.filter_struct.eta             = 1e-4; end
            if ~isfield(obj.opts.filter_struct,'lagrange_filter'), obj.opts.filter_struct.lagrange_filter = 0; end
            if ~isfield(obj.opts.filter_struct,'theta_min'),       obj.opts.filter_struct.theta_min       = 1e-4; end
            
            % maximum number of iterations
            if ~isfield(obj.opts,'max_iter'), obj.opts.max_iter                                 = 100; end

            % counter to check if we are close/almost optimal 
            % (in case of less progress for almost_opt_count iterations)
            if ~isfield(obj.opts,'almost_opt_count'), obj.opts.almost_opt_count                 = 10; end

            % method for constraint violation check
            if ~isfield(obj.opts,'con_vio_check'), obj.opts.con_vio_check                       = 'signed-distance'; end

            % user samples in case we use the sampling strategy for constraint violation checks
            if ~isfield(obj.opts,'user_sample'), obj.opts.user_sample                           = []; end

            % turn on/off display outputs
            if ~isfield(obj.opts,'verbose'), obj.opts.verbose                                   = 1; end

            % set up logger
            if ~isfield(obj.opts,'verbose') || ~obj.opts.verbose
                % no display
                obj.log = casos.package.Logger.Off;
            else
                % display debug messages
                obj.log = casos.package.Logger.Debug;
            end

            % Initialize solvers
            buildproblem(obj,nlsos);
        end

        function s = get_stats(obj)
            % Return stats.
            s = obj.sequential_stats;
            s = addfield(obj.status,s);
        end

        function s = get_info(obj)
            % Return info.
            s = get_info@casos.package.solvers.SosoptCommon(obj);
            s.sossol_convex = obj.solver_convex.get_info;
            s.sossol_constraint = obj.solver_con_vio.get_info;
            s.sossol_correction = obj.solver_soc.get_info;
        end
    end
end
