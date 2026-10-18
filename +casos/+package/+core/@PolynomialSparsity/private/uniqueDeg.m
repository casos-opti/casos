% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [coeffs,degmat,I] = uniqueDeg(coeffs,degmat)
% Make degree matrix unique and return corresponding coefficients.
%
% This function ensures that the monomials are in graded REVERSE
% lexicographic order.

nt = size(coeffs,1);
ne = size(coeffs,2);

% sort by ascending degree
degsum = sum(degmat,2);

% make degree matrix unique
[degmat2,id,ic] = unique([degsum fliplr(degmat)],'rows','sorted');

% reverse order of degrees
degmat = fliplr(degmat2(:,2:end));

if isa(coeffs,'casadi.Sparsity')
    % take union of repeated coefficients
    [ii,jj] = get_triplet(coeffs);
    coeffs = casadi.Sparsity.triplet(length(id),ne,ic(ii+1)-1,jj);

else
    % sum repeated coefficients
    summat = sparse(ic,1:nt,1,length(id),nt);
    coeffs = (summat*coeffs);
end

% undocumented: return indices
I = {id ic};

end
