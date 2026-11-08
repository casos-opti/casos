% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

classdef (Sealed) SCSInterface < casos.package.solvers.AbstractSCSInterface
% Interface for conic solver SCS.

properties (Constant, Access=protected)
    scs_options = [casos.package.solvers.ConicSolver.conic_options
        {'scs', 'Options to be passed to SCS.'}
    ];
end

methods (Static)
    function options = get_options
        % Return static options.
        options = casos.package.solvers.SCSInterface.scs_options;
    end
end

methods
    function obj = SCSInterface(name,conic,varargin)
        % Construct SCS interface.
        obj@casos.package.solvers.AbstractSCSInterface(name,conic,varargin{:});
        
        if ~isfield(obj.opts, 'scs')
            obj.opts.scs = struct;
        end


        % default options see https://www.cvxgrp.org/scs/api/settings.html
        if ~isfield(obj.opts.scs,'verbose'), obj.opts.scs.verbose = 0; end
        if ~isfield(obj.opts.scs,'use_qdldl'), obj.opts.scs.use_qdldl = 1; end
        if ~isfield(obj.opts.scs,'use_indirect'), obj.opts.scs.use_indirect = 0; end
        if ~isfield(obj.opts.scs,'dense'), obj.opts.scs.dense = 0; end
        if ~isfield(obj.opts.scs,'gpu'), obj.opts.scs.gpu = 0; end
        if ~isfield(obj.opts.scs,'normalize'), obj.opts.scs.normalize = 1; end
        if ~isfield(obj.opts.scs,'scale'), obj.opts.scs.scale = 0.1; end
        if ~isfield(obj.opts.scs,'adaptive_scale'), obj.opts.scs.adaptive_scale = 1; end
        if ~isfield(obj.opts.scs,'adaptive_diag_scale'), obj.opts.scs.adaptive_diag_scale = 1; end
        if ~isfield(obj.opts.scs,'rho_x'), obj.opts.scs.rho_x = 1e-6; end
        if ~isfield(obj.opts.scs,'max_iters'), obj.opts.scs.max_iters = 1e5; end
        if ~isfield(obj.opts.scs,'eps_abs'), obj.opts.scs.eps_abs = 1e-4; end
        if ~isfield(obj.opts.scs,'eps_rel'), obj.opts.scs.eps_rel = 1e-4; end
        if ~isfield(obj.opts.scs,'eps_infeas'), obj.opts.scs.eps_infeas = 1e-7; end
        if ~isfield(obj.opts.scs,'alpha'), obj.opts.scs.alpha = 1.5; end
        if ~isfield(obj.opts.scs,'time_limit_secs'), obj.opts.scs.time_limit_secs = 0; end
        if ~isfield(obj.opts.scs,'warm_start'), obj.opts.scs.warm_start = 0; end
        if ~isfield(obj.opts.scs,'acceleration_lookback'), obj.opts.scs.acceleration_lookback = 10; end
        if ~isfield(obj.opts.scs,'acceleration_interval'), obj.opts.scs.acceleration_interval = 5; end
        if ~isfield(obj.opts.scs,'acceleration_type_1'), obj.opts.scs.acceleration_type_1 = 1; end
        if ~isfield(obj.opts.scs,'acceleration_regularization'), obj.opts.scs.acceleration_regularization = 1e-8; end
        if ~isfield(obj.opts.scs,'acceleration_relaxation'), obj.opts.scs.acceleration_relaxation = 1.0; end
        % char* settings default to NULL, i.e. not set (empty char in MATLAB)
        if ~isfield(obj.opts.scs,'write_data_filename'), obj.opts.scs.write_data_filename = ''; end
        if ~isfield(obj.opts.scs,'log_csv_filename'), obj.opts.scs.log_csv_filename = ''; end

    end

    function s = info(obj)
        % Overwriting ConicSolver.info
        s = info@casos.package.solvers.ConicSolver(obj);
        s.scs = obj.solver_info;
    end
end

methods (Static, Access=protected)
    %% Static helper functions
    function S = sparsity_triangular(k)
        % Return lower-triangular sparsity pattern.
        S = casadi.Sparsity.lower(k);
    end
end

methods (Access=protected)
    function [x,y,s] = call_solver(obj,data,K)
        % Call SCS solver.

        % options to SCS
        opts = obj.opts.scs;

        % remove P if all entries are sparse
        % see: https://www.cvxgrp.org/scs/api/matrices.html 
        if nnz(data.P) == 0, data = rmfield(data, 'P'); end
        
        % call SCS
        [x,y,s,obj.solver_stats] = scs(data,K,opts);

        % check solver status
        status_val = obj.solver_stats.status_val;
        if status_val == -1
            % primal unbounded / dual infeasible
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,'Conic problem is dual infeasible.')
        elseif status_val == -2
            % primal infeasible / dual unbounded
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_INFEASIBLE;
            assert(~obj.opts.error_on_fail,'Conic problem is primal infeasible.')
        elseif ismember(status_val, [2 -6 -7])
            % inaccurate solution
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_NAN;
            assert(~obj.opts.error_on_fail,'Optimizer did not reach desired accuracy (Status: %s).', obj.solver_stats.status)
        elseif status_val < -2
            % failure
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_LIMITED;
            assert(~obj.opts.error_on_fail,'Optimizer failed (Status: %s).', obj.solver_stats.status)
        else
            % success
            obj.status = casos.package.UnifiedReturnStatus.SOLVER_RET_SUCCESS;
        end
    end
end

end
