% SPDX-FileCopyrightText: 2024, 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [Z,K,z,Mp,Md] = grambasis(S,I,newton_solver)
% Return Gram basis of polynomial vector.

if nargin < 3
    % no Newton simplification
    newton_solver = '';
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
I = any(Lz,1);
degmat = z.degmat(I,:);
Lz(:,~I) = [];

[Z,K,Mp,Md] = gram_internal(Lz,degmat,z.indets);	

if ~isempty(K)
    % vectorized computation of diagonal indices for each Gram block
    Kp2 = K.^2;
    % starting offset of each block
    block_offsets = [0; cumsum(Kp2(1:end-1))];   
    
    % diagonal positions within a kxk block
    diag_offsets = arrayfun(@(k) ((0:k-1)*k + (1:k)).', K, 'UniformOutput', false);
    diag_idx = vertcat(diag_offsets{:}) + repelem(block_offsets(:), K(:));
    
    % map polynomials into the Gram monomial basis Z
    poly_in_basis = poly2basis(casos.PS(S), Z);
    zero_rows     = find(full(casos.PD(poly_in_basis))==0);   
    
    % build block-diagonal "ones" pattern once
    total = sum(K);
    % build via sparse accumulation
    rows = zeros(total,1);
    cols = zeros(total,1);
    off = 0; p = 0;
    for k = K'
        r = (1:k).' + off;
        % all combinations in block: r repeated
        rows(p+1:p+k*k) = repmat(r, k, 1);
        cols(p+1:p+k*k) = repelem(r, k);
        p = p + k*k;
        off = off + k;
    end
    idx_static = sparse(rows, cols, true, total, total);
    idx = idx_static;
    
    while true
        % for each zero row, count how many nonzero entries it has in Mp
        nnz_per_row     = sum(spones(Mp*diag(idx(idx_static))), 2);
        single_nnz_rows = find(nnz_per_row==1);
        lia = ismember(single_nnz_rows, zero_rows);

        if all(lia==false); break; end

        % rows to process
        rows_to_fix = single_nnz_rows(lia);

        [~,loc] = ind2sub(size(Mp(rows_to_fix,:)), find(Mp(rows_to_fix,:)==1));

        % map back to original column indices
        loc2 = ismember(diag_idx, loc);

        % remove column (i,:) and row at (:,i)
        idx(loc2,:) = 0;
        idx(:,loc2) = 0;

        % remove monomial
        Lz(loc2) = false;
    end
    
    [Z,K,Mp,Md] = gram_internal(Lz,degmat,z.indets);
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
