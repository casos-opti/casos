% SPDX-FileCopyrightText: 2026 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

%% ------------------------------------------------------------------------
%
%
%   Short Description:  Demo for region-of-attraction estimation using 
%                       sequential quadratic SOS with a filter line-search. 
%                       Calculate an inner-estimate of the
%                       region-of-attraction for the closed-loop
%                       Van-der-Pol Oscilator.
%
%   Reference: -
%           
%
%--------------------------------------------------------------------------

% system states
x = casos.Indeterminates('x',2,1);

% system dynamics
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)];

% initial Lyapunov function
Vinit = 1.5*x(1)^2 - x(1)*x(2) + x(2)^2;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:6));

% SOS multiplier
s = casos.PS.sym('s',monomials(x,1:3),'gram');

% enforce positivity
l = 1e-6*(x'*x);
% quadratic reference function
p = x'*x/9;

%% Setup solver
% options
opts = struct; %('sossol','mosek');

% nonlinear SOS problem
sos = struct('x',[V;s],'p',[]);

% objective
sos.f = int((V-p)^2,x,[-1.5 1.5]);

% constraints
sos.g = [
            s;
            V-l; 
            s*(V-1)-nabla(V,x)*f-l;
];

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', 2);
opts.Kc = struct('sos', 3);

% build sequential solver
S = casos.nlsossol('S','sequential',sos,opts);


%% Solve nonlinear problem
sol = S('x0',[Vinit; (x'*x)]); 

%% Plot results
% Lyapunov function
Vfun = to_function(sol.x(1));
% reference function
pfun = to_function(p);

figure
fcontour(@(x,y) full(Vfun(x,y)), [-3 3], 'b', "LevelList", [1 1])
hold on
fcontour(@(x,y) full(pfun(x,y)), [-3 3], 'r--', "LevelList", [1 1])
hold off
grid on
