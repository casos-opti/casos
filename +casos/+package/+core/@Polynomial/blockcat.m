% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function p = blockcat(a,b,c,d)
% Block concatenate polynomials.

p = a.new_poly;

% check for operators
tf = [is_operator(a) is_operator(b) is_operator(c) is_operator(d)];

assert(all(~tf) || all(tf), 'Must not mix polynomials and operators.')

% concatenate coefficient matrices
[S,p.coeffs] = coeff_blkcat(a.get_sparsity, ...
                            b.get_sparsity, ...
                            c.get_sparsity, ...
                            d.get_sparsity, ...
                            a.coeffs, ...
                            b.coeffs, ...
                            c.coeffs, ...
                            d.coeffs ...
);

% set sparsity
p = set_sparsity(p,S);

end
