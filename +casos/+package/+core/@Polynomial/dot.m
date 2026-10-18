% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function c = dot(a,b)
% Scalar dot product or duality pairing.

c = a.new_poly;

if is_operator(a) && ~is_operator(b)
    % evaluate dual operator (linear form)
    assert(is_dual(a), 'Only allowed for duals. Use "evaluate" instead.')

    if ~check_sz_equal(a.sparsity_in,b)
        % dimension mismatch
        throw(casos.package.core.IncompatibleSizesError.other(a.sparsity_in,b));
    end

    % evaluation of operator on polynomial
    [S,c.coeffs] = coeff_evaluate(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);

else
    % inner product of polynomials
    assert(is_operator(a) == is_operator(b), 'Must not mix polynomials and operators.')

    if ~check_sz_equal(a,b)
        % dimension mismatch
        throw(casos.package.core.IncompatibleSizesError.other(a,b));
    end

    % compute dot product
    [S,c.coeffs] = coeff_dot(a.get_sparsity,b.get_sparsity,a.coeffs,b.coeffs);
end

% set sparsity
c = set_sparsity(c,S);

end
