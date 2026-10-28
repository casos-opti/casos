% SPDX-FileCopyrightText: 2024, 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [Z,K_red,z,Mp,Md,K] = grambasis(S,I,prune,solver)
% Return Gram basis of polynomial vector.

if nargin < 3
    % no simplification
    prune = false;
end

if prune == 2 && nargin < 4
    error('Please provide a solver to perform the Newton reduction.')
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

% build output sparsity, cone sizes and gram basis to output sparsity mapping 
switch prune
    case 0
        % No reduction algorithm
        [Z,K,Mp,Md] = gram_internal(Lz,degmat,z.indets);
        K_red = K; % no change
    case 1
        % Zero Diagonal Algorithm
        %
        % [1] Seiler, P.; Zheng, Q.; Balas, G. Simplification Methods 
        %     for Sum-of-Squares Programs. arXiv March 6, 2013.
        %
        [Z,K_red,Mp,Md,K] = zero_diagonal_alg(Lz,degmat,z.indets,S,I,idx,lp);
    case 2
        % Newton Polytope Reduction
        %
        % [2] Lofberg, J. Pre- and Post-Processing Sum-of-Squares Programs 
        %     in Practice. IEEE Trans. Automat. Contr. 2009, 54 (5), 
        %     1007–1011. https://doi.org/10.1109/TAC.2009.2017144
        %
        Ldegmat_red = Ldegmat(idx,:);
        Lz_red = arrayfun(@(i) newton_reduce(S.degmat(Ldegmat_red(i,:),Iv),degmat,solver,Lz(i,:)'), 1:lp, 'UniformOutput', false);
        Lz_red = horzcat(Lz_red{:})';
        %[Z,K,Mp,Md] = gram_internal(Lz,degmat,z.indets,Lz_red);
        [Z,K_red,Mp,Md] = gram_internal(Lz_red,degmat,z.indets);
        K = full(sum(Lz,2));
    otherwise
        error('This pruning option is not valid.')
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
