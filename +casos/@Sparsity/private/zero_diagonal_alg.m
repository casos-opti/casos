% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [Z,K_red,Mp,Md,K] = zero_diagonal_alg(Lz, degmat, indets, S, I, idx,lp)

[Z,K,Mp,Md] = gram_internal(Lz,degmat,indets);
K_red = K;

if isempty(K); return; end

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
    idx(loc2,:) = false;
    idx(:,loc2) = false;
end

% remove monomial and update (K,Z,Mp,Md)
lin = sub2ind(size(Lz), row_idx(loc2_saver), col_idx(loc2_saver));
Lz(lin) = false;

% update the mappings Mp and Md
Mp = Mp*diag(idx(idx_static));
Md = Md*diag(idx(idx_static));
% remove empty columns
empty_cols = ~any(Mp,1);
Mp(:, empty_cols) = [];
Md(:, empty_cols) = [];
% remove empty rows
empty_rows = ~any(Mp,2);
Mp(empty_rows, :) = [];
Md(empty_rows, :) = [];
% update output sparsity
[i, j, degmat] = Z.get_tuplet;
keep_mono = ~empty_rows;
Z = casos.Sparsity.tuplet(lp, 1, i(keep_mono), j(keep_mono), Z.indets, degmat(keep_mono,:));

% update the cone sizes
K_red = full(sum(Lz,2));

end
