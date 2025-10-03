% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [tf,I] = is_monom(obj)
% Check if polynomial sparsity pattern is a vector of monomial terms.

[ia,ja] = get_triplet(obj.coeffs);

% pattern is a vector of monomial terms 
% if all coefficients have exactly one nonzero entry
tf = issorted(ja,'strictascend');

% return order of monomials in vector
I = ia + 1; % NOTE: Casadi has zero-based index

end
