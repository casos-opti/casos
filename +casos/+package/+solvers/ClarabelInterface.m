% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (Sealed) ClarabelInterface < casos.package.solvers.AbstractSCSInterface
% Interface for conic solver Clarabel.

properties (Constant, Access=protected)
    clarabel_options = [casos.package.solvers.ConicSolver.conic_options
        {'clarabel', 'Options to be passed to Clarabel.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.ClarabelInterface.clarabel_options;
    end
end

methods
    function obj = ClarabelInterface(name,conic,varargin)
        % Construct Clarabel interface.
        obj@casos.package.solvers.AbstractSCSInterface(name,conic,varargin{:});

        % default options
        % see default setting structure of Clarabel MATLAB interface
           % default options
        % see default setting structure of Clarabel MATLAB interface
        if ~isfield(obj.opts, 'clarabel')
            obj.opts.clarabel = struct;
        end

        clarabel_default_opts = DefaultSettings;

        % Integer fields
        if ~isfield(obj.opts.clarabel, 'max_iter')
            obj.opts.clarabel.max_iter = clarabel_default_opts.max_iter;
        end
        if ~isfield(obj.opts.clarabel, 'equilibrate_max_iter')
            obj.opts.clarabel.equilibrate_max_iter = clarabel_default_opts.equilibrate_max_iter;
        end
        if ~isfield(obj.opts.clarabel, 'iterative_refinement_max_iter')
            obj.opts.clarabel.iterative_refinement_max_iter = clarabel_default_opts.iterative_refinement_max_iter;
        end

        % Boolean fields
        if ~isfield(obj.opts.clarabel, 'verbose')
            obj.opts.clarabel.verbose = clarabel_default_opts.verbose;
        end
        if ~isfield(obj.opts.clarabel, 'equilibrate_enable')
            obj.opts.clarabel.equilibrate_enable = clarabel_default_opts.equilibrate_enable;
        end
        if ~isfield(obj.opts.clarabel, 'direct_kkt_solver')
            obj.opts.clarabel.direct_kkt_solver = clarabel_default_opts.direct_kkt_solver;
        end
        if ~isfield(obj.opts.clarabel, 'static_regularization_enable')
            obj.opts.clarabel.static_regularization_enable = clarabel_default_opts.static_regularization_enable;
        end
        if ~isfield(obj.opts.clarabel, 'dynamic_regularization_enable')
            obj.opts.clarabel.dynamic_regularization_enable = clarabel_default_opts.dynamic_regularization_enable;
        end
        if ~isfield(obj.opts.clarabel, 'iterative_refinement_enable')
            obj.opts.clarabel.iterative_refinement_enable = clarabel_default_opts.iterative_refinement_enable;
        end
        if ~isfield(obj.opts.clarabel, 'presolve_enable')
            obj.opts.clarabel.presolve_enable = clarabel_default_opts.presolve_enable;
        end

        % Double fields
        if ~isfield(obj.opts.clarabel, 'time_limit')
            obj.opts.clarabel.time_limit = clarabel_default_opts.time_limit;
        end
        if ~isfield(obj.opts.clarabel, 'max_step_fraction')
            obj.opts.clarabel.max_step_fraction = clarabel_default_opts.max_step_fraction;
        end
        if ~isfield(obj.opts.clarabel, 'tol_gap_abs')
            obj.opts.clarabel.tol_gap_abs = clarabel_default_opts.tol_gap_abs;
        end
        if ~isfield(obj.opts.clarabel, 'tol_gap_rel')
            obj.opts.clarabel.tol_gap_rel = clarabel_default_opts.tol_gap_rel;
        end
        if ~isfield(obj.opts.clarabel, 'tol_feas')
            obj.opts.clarabel.tol_feas = clarabel_default_opts.tol_feas;
        end
        if ~isfield(obj.opts.clarabel, 'tol_infeas_abs')
            obj.opts.clarabel.tol_infeas_abs = clarabel_default_opts.tol_infeas_abs;
        end
        if ~isfield(obj.opts.clarabel, 'tol_infeas_rel')
            obj.opts.clarabel.tol_infeas_rel = clarabel_default_opts.tol_infeas_rel;
        end
        if ~isfield(obj.opts.clarabel, 'tol_ktratio')
            obj.opts.clarabel.tol_ktratio = clarabel_default_opts.tol_ktratio;
        end
        if ~isfield(obj.opts.clarabel, 'reduced_tol_gap_abs')
            obj.opts.clarabel.reduced_tol_gap_abs = clarabel_default_opts.reduced_tol_gap_abs;
        end
        if ~isfield(obj.opts.clarabel, 'reduced_tol_gap_rel')
            obj.opts.clarabel.reduced_tol_gap_rel = clarabel_default_opts.reduced_tol_gap_rel;
        end
        if ~isfield(obj.opts.clarabel, 'reduced_tol_feas')
            obj.opts.clarabel.reduced_tol_feas = clarabel_default_opts.reduced_tol_feas;
        end
        if ~isfield(obj.opts.clarabel, 'reduced_tol_infeas_abs')
            obj.opts.clarabel.reduced_tol_infeas_abs = clarabel_default_opts.reduced_tol_infeas_abs;
        end
        if ~isfield(obj.opts.clarabel, 'reduced_tol_infeas_rel')
            obj.opts.clarabel.reduced_tol_infeas_rel = clarabel_default_opts.reduced_tol_infeas_rel;
        end
        if ~isfield(obj.opts.clarabel, 'reduced_tol_ktratio')
            obj.opts.clarabel.reduced_tol_ktratio = clarabel_default_opts.reduced_tol_ktratio;
        end
        if ~isfield(obj.opts.clarabel, 'equilibrate_min_scaling')
            obj.opts.clarabel.equilibrate_min_scaling = clarabel_default_opts.equilibrate_min_scaling;
        end
        if ~isfield(obj.opts.clarabel, 'equilibrate_max_scaling')
            obj.opts.clarabel.equilibrate_max_scaling = clarabel_default_opts.equilibrate_max_scaling;
        end
        if ~isfield(obj.opts.clarabel, 'linesearch_backtrack_step')
            obj.opts.clarabel.linesearch_backtrack_step = clarabel_default_opts.linesearch_backtrack_step;
        end
        if ~isfield(obj.opts.clarabel, 'min_switch_step_length')
            obj.opts.clarabel.min_switch_step_length = clarabel_default_opts.min_switch_step_length;
        end
        if ~isfield(obj.opts.clarabel, 'min_terminate_step_length')
            obj.opts.clarabel.min_terminate_step_length = clarabel_default_opts.min_terminate_step_length;
        end
        if ~isfield(obj.opts.clarabel, 'static_regularization_constant')
            obj.opts.clarabel.static_regularization_constant = clarabel_default_opts.static_regularization_constant;
        end
        if ~isfield(obj.opts.clarabel, 'static_regularization_proportional')
            obj.opts.clarabel.static_regularization_proportional = clarabel_default_opts.static_regularization_proportional;
        end
        if ~isfield(obj.opts.clarabel, 'dynamic_regularization_eps')
            obj.opts.clarabel.dynamic_regularization_eps = clarabel_default_opts.dynamic_regularization_eps;
        end
        if ~isfield(obj.opts.clarabel, 'dynamic_regularization_delta')
            obj.opts.clarabel.dynamic_regularization_delta = clarabel_default_opts.dynamic_regularization_delta;
        end
        if ~isfield(obj.opts.clarabel, 'iterative_refinement_reltol')
            obj.opts.clarabel.iterative_refinement_reltol = clarabel_default_opts.iterative_refinement_reltol;
        end
        if ~isfield(obj.opts.clarabel, 'iterative_refinement_abstol')
            obj.opts.clarabel.iterative_refinement_abstol = clarabel_default_opts.iterative_refinement_abstol;
        end
        if ~isfield(obj.opts.clarabel, 'iterative_refinement_stop_ratio')
            obj.opts.clarabel.iterative_refinement_stop_ratio = clarabel_default_opts.iterative_refinement_stop_ratio;
        end

        % Enum fields
        if ~isfield(obj.opts.clarabel, 'direct_solve_method')
            obj.opts.clarabel.direct_solve_method = clarabel_default_opts.direct_solve_method;
        end

        if ~isfield(obj.opts.clarabel, 'chordal_decomposition_enable')
            obj.opts.clarabel.chordal_decomposition_enable = clarabel_default_opts.chordal_decomposition_enable;
        end
        if ~isfield(obj.opts.clarabel, 'chordal_decomposition_compact')
            obj.opts.clarabel.chordal_decomposition_compact = clarabel_default_opts.chordal_decomposition_compact;
        end
        if ~isfield(obj.opts.clarabel, 'chordal_decomposition_complete_dual')
            obj.opts.clarabel.chordal_decomposition_complete_dual = clarabel_default_opts.chordal_decomposition_complete_dual;
        end
        if ~isfield(obj.opts.clarabel, 'chordal_decomposition_merge_method')
            obj.opts.clarabel.chordal_decomposition_merge_method = clarabel_default_opts.chordal_decomposition_merge_method;
        end

    end

    function s = info(obj)
        % Overwriting ConicSolver.info
        s = info@casos.package.solvers.ConicSolver(obj);
        s.clarabel = obj.solver_info;
    end
end

methods (Static, Access=protected)
    %% Static helper functions
    function S = sparsity_triangular(k)
        % Return upper-triangular sparsity pattern.
        S = casadi.Sparsity.upper(k);
    end
end

methods (Access=protected)
    function [x,z,s] = call_solver(obj,data,K)
        % Call Clarabel solver.

        cones = cell(1, 1);

        % options to SCS
        opts = obj.opts.clarabel;

        % rewrite box cone into linear variables
        % (temporary workaround)
        if ~isempty(K.bl)
            % order of slacks: (sz, sl, t, sb, ...)
            A = mat2cell(data.A,[K.z+K.l 1 length(K.bl) sum(K.q)+sum(K.s*(K.s+1)/2)]);
            b = mat2cell(data.b,[K.z+K.l 1 length(K.bl) sum(K.q)+sum(K.s*(K.s+1)/2)]);

            % introduce nonnegative variables
            data.A = [A{1}; A{3};      -A{3};      A{4}];
            data.b = [b{1}; b{3}-K.bl; -b{3}-K.ub; b{4}];

            K.l = K.l + length(K.bl);
        end

        % add zero cone for the Clarabel wrapper
        if K.z > 0
            cones{end+1} = zeroConeT(K.z);
        end

        % add nonnegative cone
        if K.l > 0
            cones{end+1} = NonnegativeConeT(K.l);
        end

        % add second order cone
        if K.q > 0
            cones{end+1} = SecondOrderConeT(K.q);
        end

        % stack PSD cones
        num_psd = length(K.s);
        if num_psd > 0
            psd_cones = arrayfun(@PSDTriangleConeT, K.s, 'UniformOutput', false);
            % concatenate
            cones = [cones, psd_cones']; 
        end

        % convert cell array to a proper structure 
        cones = vertcat(cones{:});
        
        % call clarabel mex
        sol = clarabel_mex(data.P,data.c,data.A,data.b,cones,opts);
        
        % extract primal, dual and slack variables
        x = sol.x;
        z = sol.z;
        s = sol.s;
        
        % extract infos
        clarabel_info.solve_time = sol.solve_time;
        clarabel_info.status     = sol.status;
        clarabel_info.iterations = sol.iterations;
        clarabel_info.r_prim     = sol.r_prim;
        clarabel_info.r_dual     = sol.r_dual;
       
        obj.solver_stats.clarabel_info = clarabel_info;
       
        % get solution status
        if strcmp(sol.status,'Solved')
            % success
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        elseif strcmp(sol.status,'AlmostSolved')
            % success because we are feasible but did not solve until thresholds
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        elseif strcmp(sol.status,'Unsolved')
            % unsolved
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'MaxIterations')
            % max iterations reached
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'MaxTime')
            % max time reached
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'InsufficientProgress')
            % InsufficientProgress
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        elseif strcmp(sol.status,'Unknown')
            % InsufficientProgress
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_UNKNOWN;
        elseif strcmp(sol.status,'NumericalError')
            % InsufficientProgress
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
        else
            % primal unbounded / dual infeasible
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,['Conic problem is ' sol.status])
        end
    end
end

end
