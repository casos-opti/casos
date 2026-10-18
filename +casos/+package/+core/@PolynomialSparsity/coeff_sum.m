% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function [S,coeffs] = coeff_sum(obj,coeffs,dim)
% Compute sum of polynomial coefficient matrix.

sz = sizeofMatrixOp(obj,dim);

if isa(coeffs,'casadi.Sparsity')
    % compute sparsity pattern
    coeffs = sparsity_sum(coeffs,size(obj),dim);
    degmat = obj.degmat;
    indets = obj.indets;

else
    % prepare for operation on coefficients
    cfa = prepareMatrixOp(obj,coeffs,dim);

    % sum coefficients
    cfb = evalMatrixOp(cfa,@sum,dim);
    
    % finish operation on coefficients
    coeffs = finishMatrixOp(cfb,sz,dim);

    % remove zero terms
    [coeffs,degmat,indets] = removeZero(coeffs,obj.degmat,obj.indets);
end

% set dimensions of output
S = new_from_coefficients(coeffs,degmat,indets,sz);

end
