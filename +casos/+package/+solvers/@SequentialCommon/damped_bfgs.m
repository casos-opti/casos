% SPDX-FileCopyrightText: 2025 Institute of Flight Mechanics and Controls, University of Stuttgart
% SPDX-FileCopyrightText: Author(s): Torbjørn Cunis and Jan Olucak <tcunis@ifr.uni-stuttgart.de>
% SPDX-FileContributor: For a full list of contributors, see <https://github.com/ifr-ofc/casos>
%
% SPDX-License-Identifier: GPL-3.0-only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Short Description:
%
% Implementation of the damped Broyden–Fletcher–Goldfarb–Shanno (BFGS) method
% to estimate the Hessian.
% For details see:
%
% Nocedal, J., Wright, S. - Numerical Optimization, Springer, 2006, on page
% 536-537
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Bk = damped_bfgs(obj,Bk,x_k,p0,sol_iter,iter)

% update BFGS
s_val = full( obj.eval_s(x_k,sol_iter.x_k1,p0) );
y_val = full( obj.eval_y(x_k,sol_iter.x_k1,p0,sol_iter.dual_k1) );

s_dm = casadi.DM(s_val);
y_dm = casadi.DM(y_val);
B_dm = casadi.DM(Bk);

% Powell damping
sy  = full(s_dm'*y_dm);
sBs = full(s_dm'*B_dm*s_dm);

if sy >= 0.2*sBs
    theta = 1;
else
    theta = (0.8*sBs)/(sBs-sy);
end

r_val = full( obj.eval_r(Bk,theta,y_val,s_val) );

% update Bk
% B = B + (r*r')/(s'*r) - (B*s*s'*B)/(s'*B*s);
Bk =  full( obj.damped_bfgs(Bk,r_val,s_val) );

% check symmetry of Hessian and ajdust if necessary
bk_nonsym = norm(Bk-Bk',inf);

if bk_nonsym > eps
    Bk = (Bk+Bk')/2;
end

% check if Hessian is PD, if not apply perturbation
pd = true;
try
    dE = full(ldl(casadi.DM(Bk)));
    if min(dE) <= 0
        pd = false;
    end
catch
    pd = false;
end

if ~pd
    h_diag = diag(Bk);
    row_sum = sum(abs(Bk),2);
    sigma = min(h_diag - (row_sum - abs(h_diag)));
    pert = (-sigma + 2*eps*max(abs(h_diag)))/(1-2*eps);
    Bk = Bk + pert*eye(length(Bk));
end


end
