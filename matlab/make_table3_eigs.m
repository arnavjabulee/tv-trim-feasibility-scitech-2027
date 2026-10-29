function results = make_table3_eigs()
% MAKE_TABLE3_EIGS
%
% Computes Table 3 values:
%   max real eig(A_r)
%   max real eig(A_r - B_r K)
%
% Uses the same reduced state as the paper:
%   [u v w p q r phi theta psi]
%
% Uses TVC-only inputs:
%   [DELTAZ DELTAW]
%
% Assumed input order from linearize_ice_6dof:
%   [PLA DELTAZ DELTAW]

clc

db = ice_init();

% ------------------------------------------------------------
% Table 3 cases
% ------------------------------------------------------------
cases = [
    make_case("Main envelope",          6000,  250)
    make_case("Detached branch",       18000,  480)
    make_case("Near-boundary/unstable",16000,  450)
    make_case("Turn case",                 0,  300)
];

% ------------------------------------------------------------
% Reduced state:
% [u v w p q r phi theta psi]
% ------------------------------------------------------------
ix_r = [1 2 3 4 5 6 7 8 9];

% ------------------------------------------------------------
% TVC-only inputs:
% [DELTAZ DELTAW]
% ------------------------------------------------------------
iu_tvc = [2 3];

% ------------------------------------------------------------
% LQR weights
% Must match your report/code.
% This is the 9-state version, without zD.
% ------------------------------------------------------------
Q = diag([ ...
    5e0 ...  % u
    1e-1 ... % v
    1e-1 ... % w
    2e0  ... % p
    3e0  ... % q
    2e0  ... % r
    3e0  ... % phi
    3e0  ... % theta
    5e-1 ... % psi
]);

R = diag([ ...
    2e0 ... % DELTAZ
    5e0 ... % DELTAW
]);

results = struct([]);

fprintf('\nTable 3 eigenvalue results\n');
fprintf('============================================================\n');
fprintf('%-25s %8s %10s %16s %16s\n', ...
    'Case','h_m','Vt_mps','maxRe_A','maxRe_Acl');
fprintf('------------------------------------------------------------\n');

for k = 1:numel(cases)
    c = cases(k);

    try
        % Trim
        [x_trim12, u_trim, info] = trim_6dof_level_local(c.Vt, c.h, db);

        if info.exitflag <= 0 || info.resnorm > 1e-5
            warning('Trim may be poor for %s: exitflag=%g, resnorm=%.3e', ...
                c.name, info.exitflag, info.resnorm);
        end

        % Linearise full model
        [A,B] = linearize_ice_6dof(x_trim12, u_trim, db);

        % Reduced plant
        Ar = A(ix_r, ix_r);
        Br = B(ix_r, iu_tvc);

        % LQR
        K = lqr(Ar, Br, Q, R);
        Acl = Ar - Br*K;

        eig_A = eig(Ar);
        eig_Acl = eig(Acl);

        maxRe_A = max(real(eig_A));
        maxRe_Acl = max(real(eig_Acl));

        results(k).case = c.name;
        results(k).h_m = c.h;
        results(k).Vt_mps = c.Vt;
        results(k).x_trim12 = x_trim12;
        results(k).u_trim = u_trim;
        results(k).trim_info = info;
        results(k).Ar = Ar;
        results(k).Br = Br;
        results(k).K = K;
        results(k).eig_A = eig_A;
        results(k).eig_Acl = eig_Acl;
        results(k).maxRe_A = maxRe_A;
        results(k).maxRe_Acl = maxRe_Acl;

        fprintf('%-25s %8.0f %10.0f %16.6f %16.6f\n', ...
            c.name, c.h, c.Vt, maxRe_A, maxRe_Acl);

    catch ME
        warning('Failed case %s: %s', c.name, ME.message);

        results(k).case = c.name;
        results(k).h_m = c.h;
        results(k).Vt_mps = c.Vt;
        results(k).maxRe_A = NaN;
        results(k).maxRe_Acl = NaN;

        fprintf('%-25s %8.0f %10.0f %16s %16s\n', ...
            c.name, c.h, c.Vt, 'FAILED', 'FAILED');
    end
end

fprintf('============================================================\n\n');

% ------------------------------------------------------------
% LaTeX-ready rows
% ------------------------------------------------------------
fprintf('LaTeX table rows:\n');
fprintf('------------------------------------------------------------\n');

for k = 1:numel(results)
    fprintf('%s & %.0f & %.0f & %.4f & %.4f \\\\\n', ...
        results(k).case, ...
        results(k).h_m, ...
        results(k).Vt_mps, ...
        results(k).maxRe_A, ...
        results(k).maxRe_Acl);
end

fprintf('------------------------------------------------------------\n');

end

% ============================================================
% CASE HELPER
% ============================================================

function c = make_case(name, h, Vt)
c = struct();
c.name = name;
c.h = h;
c.Vt = Vt;
end

% ============================================================
% TRIM FUNCTION
% ============================================================

function [x_trim, u_trim, info] = trim_6dof_level_local(Vt_target, alt_m, db)

zD = -alt_m;

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'MaxFunctionEvaluations',3000, ...
    'MaxIterations',400);

% unknowns = [alpha_deg, theta_rad, PLA_deg, DELTAZ_deg]
x0 = [5; deg2rad(5); 90; 0];

[sol, fval, exitflag] = fsolve(@(y) trim_residuals_local(y, Vt_target, zD, db), x0, opts);

alpha_deg  = sol(1);
theta_rad  = sol(2);
PLA_deg    = sol(3);
DELTAZ_deg = sol(4);

alpha_rad = deg2rad(alpha_deg);

u = Vt_target * cos(alpha_rad);
v = 0;
w = Vt_target * sin(alpha_rad);

x_trim = [ ...
    u; ...
    v; ...
    w; ...
    0; ...
    0; ...
    0; ...
    0; ...
    theta_rad; ...
    0; ...
    0; ...
    0; ...
    zD];

u_trim = struct();
u_trim.PLA    = PLA_deg;
u_trim.DELTAZ = DELTAZ_deg;
u_trim.DELTAW = 0;

dx = ice_step_6dof(0, x_trim, u_trim, db);
gamma_actual_deg = gamma_from_state12_local(x_trim);

residual = [dx(1); dx(3); dx(5); deg2rad(gamma_actual_deg)];

info = struct();
info.alpha_deg = alpha_deg;
info.theta_deg = rad2deg(theta_rad);
info.PLA_deg = PLA_deg;
info.DELTAZ_deg = DELTAZ_deg;
info.gamma_deg = gamma_actual_deg;
info.dx_trim = dx;
info.residual = residual;
info.resnorm = norm(residual);
info.fval = fval;
info.exitflag = exitflag;

end

function R = trim_residuals_local(y, Vt_target, zD, db)

alpha_deg  = y(1);
theta_rad  = y(2);
PLA_deg    = y(3);
DELTAZ_deg = y(4);

if abs(alpha_deg) > 45 || ...
   abs(rad2deg(theta_rad)) > 60 || ...
   PLA_deg < -20 || PLA_deg > 180 || ...
   abs(DELTAZ_deg) > 60

    R = 1e6 * ones(4,1);
    return
end

alpha_rad = deg2rad(alpha_deg);

u = Vt_target * cos(alpha_rad);
v = 0;
w = Vt_target * sin(alpha_rad);

x = [ ...
    u; ...
    v; ...
    w; ...
    0; ...
    0; ...
    0; ...
    0; ...
    theta_rad; ...
    0; ...
    0; ...
    0; ...
    zD];

u_in = struct();
u_in.PLA    = PLA_deg;
u_in.DELTAZ = DELTAZ_deg;
u_in.DELTAW = 0;

dx = ice_step_6dof(0, x, u_in, db);
gamma_actual_rad = deg2rad(gamma_from_state12_local(x));

R = [ ...
    dx(1); ...
    dx(3); ...
    dx(5); ...
    gamma_actual_rad];

if any(~isfinite(R))
    R = 1e6 * ones(4,1);
end

end

% ============================================================
% FLIGHT PATH ANGLE
% ============================================================

function gamma_deg = gamma_from_state12_local(x12)

u = x12(1);
v = x12(2);
w = x12(3);

phi   = x12(7);
theta = x12(8);
psi   = x12(9);

cphi = cos(phi);   sphi = sin(phi);
cth  = cos(theta); sth  = sin(theta);
cpsi = cos(psi);   spsi = sin(psi);

R_b2n = [ ...
    cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
    cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
    -sth,      sphi*cth,                   cphi*cth];

Vn = R_b2n * [u; v; w];

gamma_deg = rad2deg(atan2(-Vn(3), sqrt(Vn(1)^2 + Vn(2)^2)));

end