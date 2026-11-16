% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Renato Loureiro <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

% Demonstrate monomial basis reduction through zero diagonal algorithm

x = casos.PS('x', 3, 1);
f = casos.PS.sym('f');

g = [12+x(2)^2-2*x(1)^3*x(2)+2*x(2)*x(3)^2+x(1)^6-2*x(1)^3*x(3)^2+x(3)^4+x(1)^2*x(2)^2+f;
    (x(1)-1)^2+x(1)^4+f-2;
    3*x(1)^4-2*x(1)^2*x(2)+7*x(1)^2-4*x(1)*x(2)+4*x(2)^2+f;
    x(1)^4 + x(2)^4 - 2*x(1)^2*x(2)^2];

sos = struct('x', f, 'g', g, 'f', f);

% states + constraint are SOS cones
opts.Kx.lin = 1;
opts.Kx.sos = 0;
opts.Kc.sos = 4;

% ignore infeasibility
opts.error_on_fail = false;

% options in monomial pruning:
%   "none"             - no pruning
%   "zero-diagonal"    - make use the zero diagonal algorithm for reduction 
%   "newton-polytope"  - use newton polytope to reduce required monomials (default)
opts.grambasis_reduction = 'zero-diagonal';

% build the solver
S = casos.sossol('S','mosek',sos,opts); % build problem

% run solver
sol = S();
fprintf('sol.f = %d \n', full(sol.f)); % the value should be 1.7107

% obtain an accurate computation times
tbld = timeit(@() casos.sossol('S','mosek',sos,opts));
tsol = timeit(@() S());

fprintf('build time: %ds \n', tbld);
fprintf('solve time: %ds \n', tsol);
fprintf('reduction: %s \n', opts.grambasis_reduction);
fprintf('Gram basis sizes:\n')
disp(S.info.gram.Kc')
