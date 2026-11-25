% SPDX-FileCopyrightText: 2024 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

% Projection onto sum-of-squares cone

rng(0)

% indeterminate variable
x = casos.Indeterminates('x');
% some random polynomial
p = casos.PD(monomials(x,0:4),randn(5,1));

% Gram decision variable
s = casos.PS.sym('q',grambasis(p));

% projection error
e = s - p;

% define Q-SOS problem:
%   min ||s-p||^2 s.t. s is SOS
sos = struct('x',s,'f',dot(e,e));
% states is scalar SOS cone
opts = struct('Kx',struct('sos',1));

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);
% evaluate
sol = S();

fprintf('Distance to SOS cone is %g.\n', full(sol.f))
