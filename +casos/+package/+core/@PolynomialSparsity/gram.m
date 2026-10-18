% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [Z,K,Mp,Md] = gram(w)
% Return a Gram basis for a given monomial sparsity pattern.

[degmat,Lw] = get_degmat(w);

[Z,K,Mp,Md] = gram_internal(Lw,degmat,w.indets);

end
