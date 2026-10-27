% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

function buildproblem(obj,nlsos)
% Build convex subprobems for sequential algorithms.
t_build = tic;

opts = obj.opts;

% sparsity patterns
base_x = sparsity(nlsos.x);

% Jacobians
Df = jacobian(nlsos.f,nlsos.x);
Dg = jacobian(nlsos.g,nlsos.x);

% Hessian /BFGS matrix
B = casadi.SX.sym('B',nnz(base_x),nnz(base_x));

Hf = casos.PSOperator(B,base_x,base_x);

%% build quadratic SOS problem
sos_qp = [];
% convex decision variables
sos_qp.x = casos.PS.sym('x_convex',base_x);
% search direction
d = sos_qp.x - nlsos.x;
% quadratic cost function
sos_qp.f = (1/2)*dot2(Hf,d,d) + dot(Df,d);

% linear constraints
sos_qp.g = nlsos.g + dot(Dg,d);

% parameters to subproblem
% * parameters to nonlinear problem, current (nonlinear) solution, Hessian
sos_qp.p = [nlsos.p; nlsos.x; B(:)];

% SOS options
sosopt               = opts.sossol_options;
sosopt.Kx            = opts.Kx;
sosopt.Kc            = opts.Kc;
sosopt.error_on_fail = false;

sos_qp.derivatives.Hf = Hf;
sos_qp.derivatives.Df = Df + dot(casos.PSOperator.from_primal(d),Hf); % correct for shift in variables
sos_qp.derivatives.Dg = Dg;

% initialize convex SOS solver
obj.solver_convex = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_qp,sosopt);

% store basis
obj.sparsity_x  = obj.solver_convex.sparsity_x;
obj.sparsity_xl = obj.solver_convex.sparsity_xl;
obj.sparsity_xs = obj.solver_convex.sparsity_xs;
obj.sparsity_p  = sparsity(nlsos.p);
obj.sparsity_f  = sparsity(nlsos.f);
obj.sparsity_g  = obj.solver_convex.sparsity_g;
obj.sparsity_gl = obj.solver_convex.sparsity_gl;
obj.sparsity_gs = obj.solver_convex.sparsity_gs;

% basis of decision variables and parameter (needed in either branch below)
basis_nlsos_x = poly2basis(nlsos.x);
basis_nlsos_p = poly2basis(nlsos.p);

if strcmpi(obj.opts.con_vio_check,'signed-distance')
    %% constraint violation check via signed distance
    % sparsity patterns of nonlinear constraints

    % get cone dimensions
    Ml = get_dimension(obj.get_cones,opts.Kc,'lin');
    Ms = get_dimension(obj.get_cones,opts.Kc,'sos');

    % create logical indices for SOS
    I = [false(Ml,1); true(Ms,1)];

    % create function to check the linear equality contraints violation
    coordinates_lin = poly2basis(nlsos.g,~I);
    nnz_lin = coordinates_lin.nnz;

    % set upper and lower bounds
    ub = casadi.SX.sym('ub', nnz_lin, 1);
    lb = casadi.SX.sym('lb', nnz_lin, 1);

    % create function to evaluate the violation of the linear cone 
    obj.linvio = casadi.Function('linvio', {basis_nlsos_x, basis_nlsos_p, lb, ub}, {[coordinates_lin-ub; lb-coordinates_lin]});

    % get gram half-basis for nonlinear constraints
    [~,~,z] = grambasis(nlsos.g,I);

    % build unit vectors
    base_s0 = gramunit(z);

    r  = casos.PS.sym('r',sum(I));
    s0 = casos.PD(base_s0);

    sos_con_vio = [];
    % decision variables
    sos_con_vio.x = r;
    % cost function
    sos_con_vio.f = sum(r);
    % constraints
    sos_con_vio.g = nlsos.g(I) + r.*s0;
    % parameters to subproblem
    % * parameters to nonlinear problem, current (nonlinear) solution
    sos_con_vio.p = [nlsos.p; nlsos.x];

    % SOS options
    sosopt               = opts.sossol_options;
    sosopt.Kx            = struct('lin',length(r));
    sosopt.Kc            = opts.Kc;
    if isfield(sosopt.Kc, 'lin')
        sosopt.Kc = rmfield(sosopt.Kc, 'lin');
    end
    sosopt.error_on_fail = true;

    % initialize convex SOS solver
    obj.solver_con_vio = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_con_vio,sosopt);
else
    % evaluate constraints on sample points
    X  = casadi.SX.sym('x',nlsos.g.nvars,1);
    G = subs(nlsos.g, nlsos.g.indeterminates, X);

    obj.eval_constraint_samples = casadi.Function('g', {basis_nlsos_x, X}, {casadi.SX(G)});
end

%% Second-order correction (soc)
if obj.opts.enable_soc
    % parameter of current qp solution
    x_star_par   = casos.PS.sym('x_star',base_x);
    % constraint function to be evaluated at xk+d
    con_fun       = casos.Function('f',{nlsos.x, nlsos.p},{nlsos.g});

    sos_soc = [];
    % convex decision variables
    sos_soc.x = casos.PS.sym('x_soc',base_x);
    % search direction
    d = sos_soc.x - nlsos.x;
    % quadratic cost function
    sos_soc.f = (1/2)*dot2(Hf,d,d) + dot(Df,d);
    % augmented linear constraints (see Nocedal p.544)
    dk        = x_star_par - nlsos.x ;
    sos_soc.g = dot(Dg,d) + con_fun(x_star_par,nlsos.p) - dot(Dg, dk);
    % parameters to subproblem
    % * parameters to nonlinear problem, current (nonlinear) solution, Hessian,
    % solution QP
    sos_soc.p = [nlsos.p; nlsos.x; B(:);x_star_par];

    sos_soc.derivatives.Hf = Hf;
    sos_soc.derivatives.Df = Df + dot(casos.PSOperator.from_primal(d),Hf);
    sos_soc.derivatives.Dg = Dg;

    % SOS options
    sosopt    = opts.sossol_options;
    sosopt.Kx = opts.Kx;
    sosopt.Kc = opts.Kc;
    sosopt.error_on_fail = false;

    % initialize adapted convex SOS solver for soc
    obj.solver_soc = casos.package.solvers.sossolInternal('SOS',opts.sossol,sos_soc,sosopt);
end

%% Damped BFGS (see Nocedal p.536/537)
lam_gs    = casos.PS.sym('lam_gs', obj.sparsity_g);
x_k1      = casos.PS.sym('x_k1',base_x);

lagrangian = nlsos.f + dot(lam_gs,nlsos.g);

basis_xk1    = poly2basis(x_k1);
basis_lam_gs = poly2basis(lam_gs);

% Langrangian casos function for evaluation
obj.eval_L = casos.Function('L',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, {lagrangian });

% Hessian of Langrangian; in case we use BFGS, we do not have to explictly
% compute the Hessian. The hessian call also returns the gradient of the
% Lagrangian, hence it is projected and differentiated only once.
if ~strcmpi(obj.opts.hessian_approx,'BFGS')
    % current work around, because poly. hessian does not output jacobian
    [L_xx, dLdx0] = hessian(poly2basis(lagrangian), basis_nlsos_x);
    obj.eval_Hessian = casos.Function('H',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, { L_xx });
else
    dLdx0 = jacobian(poly2basis(lagrangian), basis_nlsos_x)';
end
% search direction
obj.eval_s = casos.Function('s',{basis_nlsos_x,basis_xk1,basis_nlsos_p}, { poly2basis( x_k1 - nlsos.x)});

% delta gradient of Langrangian
dLdx = casos.Function('dLdx',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, {dLdx0 });
obj.eval_dLdx = dLdx;

% y =  \nabla_x L (x_k1,lambda_k1) - \nabla_x L (x_k,lambda_k1)
obj.eval_y = casos.Function('y',{basis_nlsos_x,basis_xk1,basis_nlsos_p,basis_lam_gs}, {dLdx(basis_xk1,basis_nlsos_p,basis_lam_gs) - dLdx(basis_nlsos_x,basis_nlsos_p,basis_lam_gs) });

% using MX instead of SX is faster to build; no significant difference in
% execution time during online execution (but just a few examples)
theta = casadi.MX.sym('theta');
s     = casadi.MX.sym('s',[size(B,1),1]);
y     = casadi.MX.sym('y',[size(B,1),1]);
r     = casadi.MX.sym('r',[size(B,1),1]);
B     = casadi.MX.sym('B',[size(B,1),size(B,1)]);

% Powell-damping
obj.eval_r = casos.Function('f',{B,theta,y,s}, {theta*y+(1-theta)*B*s});

% B = B + (r*r')/(s'*r) - (B*s*s'*B)/(s'*B*s)
obj.damped_bfgs =  casadi.Function('f',{B,r,s}, {B - (B*(s*s')*B)/(s'*B*s) + (r*r')/(s'*r)});

%% evaluation functions
% nonlinear cost function
obj.eval_cost      = casos.Function('f',{basis_nlsos_x,basis_nlsos_p}, { full(nlsos.f) });
% gradient cost
obj.eval_gradCost  = casos.Function('f',{basis_nlsos_x,basis_nlsos_p}, { op2basis( Df ) });
% norm gradient langrangian
obj.eval_gradLag =  casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, { norm(dLdx0,inf ) });
obj.eval_gradLag_inf_norm  =  casos.Function('f',{basis_nlsos_x,basis_nlsos_p,basis_lam_gs}, { norm(dot(lam_gs,nlsos.g),inf) });

% step norms of primal and dual iterates
xk_step  = casadi.SX.sym('x_k',nnz(base_x),1);
xk1_step = casadi.SX.sym('x_k1',nnz(base_x),1);
dk_step  = casadi.SX.sym('dual_k',obj.sparsity_g.nnz,1);
dk1_step = casadi.SX.sym('dual_k1',obj.sparsity_g.nnz,1);
obj.eval_step_norms = casadi.Function('f',{xk_step,xk1_step,dk_step,dk1_step}, ...
    {norm(xk_step,inf), norm(xk1_step-xk_step,inf), norm(dk_step,inf), norm(dk1_step-dk_step,inf)});

% thresholds of the first-order optimality termination checks
f_val   = casadi.SX.sym('f');
gln_val = casadi.SX.sym('gln');
gl_val  = casadi.SX.sym('gl');
dxi_val = casadi.SX.sym('dxi');
tol_val = casadi.SX.sym('tol');
obj.eval_termination = casadi.Function('f',{f_val,gln_val,gl_val,dxi_val,tol_val}, ...
    {(tol_val*fmax(1,abs(f_val)) + gln_val)/dxi_val, (10*tol_val*fmax(1,abs(f_val)) + gln_val)/dxi_val, tol_val*fmax(1,gl_val)});

% directional derivative of the cost along the search direction
x_star_dir = casos.PS.sym('x_star',base_x);
obj.eval_nabla_f_dir = casos.Function('f',{basis_nlsos_x, poly2basis(x_star_dir), basis_nlsos_p}, { dot(Df, x_star_dir - nlsos.x) });

obj.eval_constraint = casos.Function('f',{basis_nlsos_x,basis_nlsos_p}, {poly2basis(nlsos.g)});

%% get parameter that are needed for initilization
obj.init_para.no_dual_var   = obj.sparsity_g.nnz;
obj.init_para.size_B        = nnz(base_x);

if strcmpi(obj.opts.con_vio_check,'signed-distance')
    obj.init_para.con_vio.no_con = length(sos_con_vio.g);
end

%% get params for display output (problem size etc)
obj.display_para.no_decVar         = obj.solver_convex.sparsity_x.nnz;
obj.display_para.no_sosCon         = length(nlsos.g);
obj.display_para.solver            = obj.opts.sossol;
obj.display_para.solver_build_time = toc(t_build);

end
