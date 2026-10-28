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

% To disable simplification:
%   - Set opts.prune to 0
opts.prune = 1; 

% Build the solver
tic
S = casos.sossol('S','mosek',sos,opts); % build problem
fprintf('build time: %ds \n', toc);

% Run solver
sol = S();
fprintf('sol.f = %d \n', full(sol.f)); % the value should be 1.7107

% Run the solver to obtain an accurate computation time
tsol = timeit(@()S());
fprintf('reduction = %d \n', opts.prune);
fprintf('time: %ds \n',tsol);
disp(S.info.gram.Kc')