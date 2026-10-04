% SPDX-FileCopyrightText: 2024, 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [Z,K,z,Mp,Md,z_rem] = grambasis(S,I,prune)
% Return Gram basis of polynomial vector.

if nargin < 3
    % no simplification
    prune = false;
end

if nargin < 2 || isempty(I)
    lp = numel(S);
    I = true(lp,1);
    idx = 1:lp;
else
    lp = nnz(I);
    idx = find(I);
end

% get logical maps for degrees and indeterminate variables
% -> Ldeg(i,j) is true iff p(i) has terms of degree(j)
% -> Lvar(i,j) is true iff p(i) has terms in indets(j)
[degree,Ldeg] = get_degree(S,I);
[indets,Lvar] = get_indets(S,I);
% remove non-indexed logicals
Ldeg(~I,:) = []; Lvar(~I,:) = [];
% detect degrees and variables
Id = any(Ldeg,1); Iv = any(Lvar,1);
% remove unused degrees or variables
Ldeg(:,~Id) = []; Lvar(:,~Iv) = [];

% min and max degree of Gram basis vector
mndg = floor(min(degree)/2);
mxdg =  ceil(max(degree)/2);

% ensure min and max degree are even
[degree,ii] = unique([degree,(2*mndg):(2*mxdg)]);
ld = length(degree); tmp = [Ldeg false(lp,ld)];
Ldeg = tmp(:,ii);

% split into even and odd monomials
Ldeg_e = Ldeg(:,mod(degree,2)==0);
Ldeg_o = Ldeg(:,mod(degree,2) >0);
% add even degrees if necessary
Ldeg_e(:,1:end-1) = Ldeg_e(:,1:end-1) | Ldeg_o;
Ldeg_e(:,2:end) = Ldeg_e(:,2:end) | Ldeg_o;

% get half-degree vector of monomials
% TODO: compute degree matrix directly (avoid unnecessary checks)
z = to_vector(casos.Sparsity.scalar(indets,mndg:mxdg));
% compute logical map for half-degree monomials
% -> Lz_deg(i,j) is true iff z(i) is of degree(j)/2
% -> Lz_var(i,j) is true iff z(i) includes indets(j)
[~,Lz_deg] = get_degree(z);
[~,Lz_var] = get_indets(z);

% build logical map for half-degree monomials
% -> Lz(i,j) is true iff Gram form of p(i) includes z(j)
% that is, z(j) is of a degree in the Gram basis for p(i)
% AND z(j) only has indeterminate variables that p(i) has, too
Lz = (Ldeg_e * Lz_deg' & ~(~Lvar * Lz_var'));

% discard monomials based on simple checks
[~,Ldegmat] = get_degmat(S);
lz = numel(z);
% perform checks vector-wise
MX = arrayfun(@(i) repmat(max(S.degmat(Ldegmat(i,:),Iv),[],1),lz,1), idx, 'UniformOutput', false);
MN = arrayfun(@(i) repmat(min(S.degmat(Ldegmat(i,:),Iv),[],1),lz,1), idx, 'UniformOutput', false);
% Imx = any(  ceil(max(p.degmat,[],1)/2) < z.degmat , 2 );
% Imn = any( floor(min(p.degmat,[],1)/2) > z.degmat , 2 );
zdm = repmat(z.degmat,lp,1);
Irem = [ceil(vertcat(MX{:})/2) < zdm, floor(vertcat(MN{:})/2) > zdm];
% Lz(:,Imx | Imn) = false;
Lz(reshape(any(Irem,2),lz,lp)') = false;

% remove unused monomials from base vector
Ir = any(Lz,1);
degmat = z.degmat(Ir,:);
Lz(:,~Ir) = [];

[Z,K,Mp,Md] = gram_internal(Lz,degmat,z.indets);	

% empty sparsity for removed monomials
z_rem = casos.Sparsity;

% apply zero diagonal algorithm
if ~isempty(K) && prune
    % indexes for later update of monomial basis
    [col_idx,row_idx] = find(Lz'==1);

    % diagonal positions within a kxk block
    diag_offsets = arrayfun(@(k) ((0:k-1)*k + (1:k)).', K, 'UniformOutput', false);
    diag_idx = vertcat(diag_offsets{:}) + reshape(repelem([0; cumsum(K(1:end-1).^2)], K(:)), [], 1);
    
    [Smat, SLmat] = get_degmat(S,I);
    SLmat = SLmat(:,any(SLmat(I,:),1));
    [Zmat, ZLmat] = get_degmat(Z);

    poly_in_basis = arrayfun(@(i) ismember(Zmat(ZLmat(i,:),:), Smat(SLmat(i,:),:), 'rows'), 1:size(idx,1), 'UniformOutput', false);
    zero_rows = find(vertcat(poly_in_basis{:})==0);

    % block-diagonal of ones, one k×k dense block per element of K
    c = repelem(1:numel(K), K);          % block id per row/col
    idx = sparse((c == c.'));            % block-diagonal mask
    idx_static = idx;

    loc2_saver = false(size(diag_idx));
    for iter = 1:sum(K)
        % for each zero row, count how many nonzero entries it has in Mp
        Mp_red = Mp*diag(idx(idx_static));
        rows_to_fix = intersect(find(sum(Mp_red,2)==1), zero_rows);

        % in case no zero equality if found on a diagonal element
        if isempty(rows_to_fix); break; end

        % map back to original column indices
        [~,loc] = find(Mp_red(rows_to_fix,:));
        loc2 = ismember(diag_idx, loc);
        loc2_saver = loc2 | loc2_saver;

        % remove column (i,:) and row at (:,i)
        idx(loc2,:) = false;    idx(:,loc2) = false;
    end

    % save removed monomials from basis (for debug)
    if nargout == 6
        Lz_del = sparse(row_idx(loc2_saver), col_idx(loc2_saver), 1, size(Lz,1), size(Lz,2));
        [i,j] = find(Lz_del');
        coeffs = casadi.Sparsity.triplet(size(Lz_del,2),lp,i-1,j-1);
        z_rem = casos.Sparsity;
        [z_rem.coeffs,z_rem.degmat] = uniqueDeg(coeffs,degmat);
        z_rem.indets = indets;
        z_rem.matdim = [lp 1];
    end

    % remove monomial and update (K,Z,Mp,Md)
    lin = sub2ind(size(Lz), row_idx(loc2_saver), col_idx(loc2_saver));
    Lz(lin) = false;
    
    % update the mappings Mp and Md 
    Mp = Mp*diag(idx(idx_static));
    Md = Md*diag(idx(idx_static));
    spar = ~any(Mp,1);
    Mp(:, spar) = [];
    Md(:, spar) = [];

    % update the cone sizes
    K = full(sum(Lz,2));
end

% build half-basis for each element
[i,j] = find(Lz');
coeffs = casadi.Sparsity.triplet(size(Lz,2),lp,i-1,j-1);
% set output
z = casos.Sparsity;
[z.coeffs,z.degmat] = uniqueDeg(coeffs,degmat);
z.indets = indets;
z.matdim = [lp 1];

end
