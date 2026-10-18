% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function Lz = newton_reduce(Pdegmat,Zdegmat,solver,Lz)
% Removes monomials outside half Newton polytope.
%
%   Inputs:
%       Pdegmat - Degree matrix of polynomial basis (each row is a monomial)
%       Zdegmat - Degree matrix of monomials to test (each row is a monomial)
%       solver  - Solver for linear programming (e.g., 'sedumi', 'mosek')
%       Lz      - Initial logical vector
%
%   Output:
%       Lz      - Sparse logical vector indicating kept monomials
%
%   Strategy inspired by: "Simplification Methods for Sum-of-Squares Programs",
%   Peter Seiler et al.

% tolerance for feasibility
tolerance = 1e-6;

% get dimensions
[nP, ~] = size(Pdegmat);
nZ = size(Zdegmat, 1);

% pre-allocate keep array (start with provided Lz)
keep = Lz;

% check if any monomial in Zdegmat*2 matches a row in Pdegmat
Z2 = Zdegmat * 2;
keep_trivial = ismember(Z2, Pdegmat, 'rows');

% early return if all are trivially kept
if all(keep_trivial)
    Lz = sparse(keep);
    return;
end

% fixed part of LP constraints (common for all iterations)
bfixed = [zeros(nP, 1); 1];

% pre-allocate solver structure (will be reused)
Slin = [];

% check remaining monomials via linear programming
for i = 1:nZ
    % skip if already kept or trivially not kept
    if keep_trivial(i) || ~keep(i)
        continue;
    end
        
    % current monomial to test (scaled by 2)
    q = Z2(i, :);
        
    % extract A matrix and b vector
    A = [Pdegmat, -ones(nP, 1); q, -1]; % linear constraint matrix
    b = bfixed;                         % right-hand side
        
    % construct augmented matrix for LP
    % variables: [x; t; s] where:
    %   x - coefficients for hyperplane
    %   t - scalar variable
    %   s - slack variables for inequalities
    n_x = size(A, 2);  % number of x variables
    n_s = size(A, 1);  % number of slack variables
        
    % build LP in form: min c'*[x; t; s] subject to a*[x; t; s] <= [b; 0]
    a = [A, -eye(n_s), zeros(n_s,1);%;        % A*[x;t] <= b
         -q, 1, zeros(1, n_s), -1];           % -q*x + t <= 0
    c = [zeros(n_x, 1); ones(n_s, 1); 1];     % objective: sum(s)
        
    % reuse or create solver instance
    if isempty(Slin) || any(size(lp_struct.a) ~= size(a))
        % create new LP structure
        lp_struct = struct('g', casadi.Sparsity.dense(size(c)), ...
                           'a', casadi.Sparsity.dense(size(a)));
        Slin = casos.conic('S', solver, lp_struct);
    end

    % Solve LP
    try
        sol = Slin('h',   sparse(size(a, 2), size(a, 2)),   ...
                   'g',   c,                                ...
                   'a',   a,                                ...
                   'lba', -inf,                             ...
                   'uba', [b; 0],                           ...
                   'lbx', [-inf(n_x, 1); zeros(n_s+1, 1)],  ...
                   'ubx', inf,                              ...
                   'x0',  sparse(size(a, 2), 1));
    catch ME
        % Solver failed
        error('Solver failed during newton Polytope reduction: %s', ME.message)
    end

    % check solver status
    status = Slin.stats.UNIFIED_RETURN_STATUS;
    if ~isequal(status,'SOLVER_RET_SUCCESS')
        continue;
    end
    
    % check feasibility (cost < tolerance means feasible)
    cost = full(sol.cost);
    if cost < tolerance
        % extract solution and identify removed monomials
        x_sol = full(sol.x(1:n_x-1));
                
        % find monomials that violate the hyperplane inequality
        % 2*z*a <= b1  ->  check if 2*z*a - x_sol(end) > 0
        b1 = full(sol.x(n_x));
        violation = Z2 * x_sol - b1 > 0;
                
        % remove violating monomials
        keep(violation) = false;
    end

end

% return as sparse logical vector
Lz = sparse(keep);

end
