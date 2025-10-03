% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,I] = restrict_terms(S,deg)
% Restrict monomials terms to specified degrees.

% choose terms to keep
I = ismember(S.degsum,deg);

% restrict coefficients
coeffs = sub(S.coeffs,find(I)-1,(1:S.numel)-1);
% restrict degree matrix
degmat = S.degmat(I,:);

% remove zero terms
[coeffs,degmat,indets] = removeZero(coeffs,degmat,S.indets);

% new sparsity pattern
S = new_from_coefficients(coeffs,degmat,indets,size(S));

end
