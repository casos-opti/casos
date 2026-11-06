% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [Z,K_red,Mp,Md,K] = newton_reduce(Lz, degmat, indets, S, I, idx, lp, Iv, solver)
% Removes monomials outside half Newton polytope.
%
%   Strategy inspired by: 
%   [1] "Simplification Methods for Sum-of-Squares Programs", Peter Seiler et al.
%   [2] "Pre- and Post-Processing Sum-of-Squares  Programs in Practice", Johan Löfberg

% get the original (non-pruned) psd sizes
K = full(sum(Lz,2));

if isempty(K)
    [Z,K_red,Mp,Md] = gram_internal(Lz,degmat,indets);
    return; 
end

% tolerance for feasibility
tolerance = 1e-6;

% get the degree matrix of the full monomial basis S 
[~,Ldegmat] = get_degmat(S);
Ldegmat_red = Ldegmat(idx,:);

Zdegmat = degmat;

% pre-compute per-j quantities
Lz_cell = cell(1,lp);           % will hold the reduced (kept) monomials per poly
Pdegmat_cell = cell(1,lp);      % monomials appearing in polynomial j (in variables Iv)
Z2_cell = cell(1,lp);           % doubled Z degree matrix (Gram monomials)
keep_cell = cell(1,lp);         % boolean mask of monomials kept for poly j
keep_trivial_cell = cell(1,lp); % monomials trivially kept (present in P)
bfixed_cell = cell(1,lp);       % RHS vector for the LP constraint
skip_j = false(1,lp);           % flag: poly j needs no reduction

nZ = size(Zdegmat, 1);

for j = 1:lp
    % monomials of polynomial j, projected onto variables Iv
    Pdegmat_cell{j} = S.degmat(Ldegmat_red(j,:),Iv);
    Lz_j = Lz(j,:)';
    
    [nP, ~] = size(Pdegmat_cell{j});
    
    % initialize the keep-mask from the current Lz row
    keep_cell{j} = Lz_j;
    
    % doubled degree matrix
    Z2_cell{j} = Zdegmat * 2;
    keep_trivial_cell{j} = ismember(Z2_cell{j}, Pdegmat_cell{j}, 'rows');
    
    % if every Gram monomial is trivially needed, skip this polynomial
    if all(keep_trivial_cell{j})
        skip_j(j) = true;
    else
        % constraint RHS: [0; ...; 0; 1] (last entry = normalization)
        bfixed_cell{j} = [zeros(nP, 1); 1];
    end
end

% single flattened loop over (j, i) pairs
Slin = [];

% build the full list of (j,i) pairs to process
[jj, ii] = ndgrid(1:lp, 1:nZ);
jj = jj(:);
ii = ii(:);

% check remaining monomials via linear programming
for k = 1:numel(jj)
    j = jj(k);
    i = ii(k);

    % Nothing to do for this poly or this monomial
    if skip_j(j)
        continue;
    end

    keep = keep_cell{j};
    keep_trivial = keep_trivial_cell{j};

    % skip if already kept or trivially not kept
    if keep_trivial(i) || ~keep(i)
        continue;
    end

    Pdegmat = Pdegmat_cell{j};
    Z2 = Z2_cell{j}; 
    bfixed = bfixed_cell{j};

    [nP, ~] = size(Pdegmat);
    q = Z2(i, :);   % doubled degree of candidate monomial i
        
    % LP feasibility problem: does there exist a nonnegative combination
    % of P's monomials plus q equal to the normalized constraint?
    % A = [Pdegmat, -1; q, -1], b = [0;...;0; 1]
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
        keep_cell{j} = keep;
    end

end

% collect the reduced monomial sets per polynomial into a single matrix
for j = 1:lp
    Lz_cell{j} = sparse(keep_cell{j});
end
Lz_red = horzcat(Lz_cell{:})';

% rebuild the Gram matrix factorization for the reduced basis
[Z,K_red,Mp,Md] = gram_internal(Lz,degmat,indets,Lz_red);

% post-process output sparsity: remove rows of Mp/Md that became empty,
% unless the corresponding monomial is still part of the original SOS basis
empty_rows  = ~any(Mp,2);
[Smat, SLmat] = get_degmat(S,I);
SLmat = SLmat(:,any(SLmat(I,:),1));
[Zmat, ZLmat] = get_degmat(Z);
poly_in_basis_cell = arrayfun(@(i) ismember(Zmat(ZLmat(i,:),:), Smat(SLmat(i,:),:), 'rows'), 1:size(idx,1), 'UniformOutput', false);
poly_in_basis = vertcat(poly_in_basis_cell{:});
remove_rows = full(empty_rows) & (~poly_in_basis);
Mp(remove_rows, :) = [];
Md(remove_rows, :) = [];

% rebuild the monomial sparsity object Z with only the kept monomials
[i, j, degmat] = Z.get_tuplet;
keep_mono = ~remove_rows;
Z = casos.Sparsity.tuplet(lp, 1, i(keep_mono), j(keep_mono), Z.indets, degmat(keep_mono,:));

end
