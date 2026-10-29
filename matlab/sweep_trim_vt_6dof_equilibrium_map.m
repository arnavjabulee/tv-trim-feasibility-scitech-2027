function results = sweep_trim_vt_6dof_equilibrium_map(mode)
% SWEEP_TRIM_VT_6DOF_EQUILIBRIUM_MAP
%
% Computes a 6-DOF level-flight equilibrium map from scratch using:
%   - normal ice_step_6dof
%   - PLA as a trim variable
%   - DELTAZ as a trim variable
%   - engine behaviour supplied by ice_block_engine
%
% This version is intended for use with an ice_block_engine that has been
% modified to extrapolate thrust/drag in PLA beyond the original table max.
%
% Unknowns:
%   y = [alpha_deg; theta_rad; PLA_deg; DELTAZ_deg]
%
% Residuals:
%   dx(1)
%   dx(3)
%   dx(5)
%   gamma_actual - gamma_target
%
% Modes:
%   'level' only

if nargin < 1
    mode = 'level';
end

mode = lower(string(mode));

db = ice_init();

% ---------------- USER OPTIONS ----------------
Vlist = 30:2:36;              % m/s
alt_list_m = 0:100:20000;      % m

showConsoleProgress = true;
trim_tol = 1e-5;
saveCSV = true;
makePlots = true;

% Allow solver to search above original engine table max if engine
% extrapolation has been enabled inside ice_block_engine.
PLA_min = min(db.eng.PLA);
PLA_max = max(db.eng.PLA) + 1000;   % adjust if needed

% Broad search bound only, not a feasibility limit from the study
DELTAZ_abs_bound_deg = 50;

nConsecutiveNoEqStop = 1;
% ---------------------------------------------

results = struct();
results.settings = struct( ...
    'Vlist', Vlist, ...
    'alt_list_m', alt_list_m, ...
    'trim_tol', trim_tol, ...
    'PLA_min', PLA_min, ...
    'PLA_max', PLA_max, ...
    'DELTAZ_abs_bound_deg', DELTAZ_abs_bound_deg, ...
    'nConsecutiveNoEqStop', nConsecutiveNoEqStop);

if mode ~= "level"
    error('This version currently supports only mode = "level".');
end

results.level = run_level_case_6dof_equilibrium_map( ...
    db, Vlist, alt_list_m, ...
    PLA_min, PLA_max, DELTAZ_abs_bound_deg, ...
    showConsoleProgress, trim_tol, nConsecutiveNoEqStop);

if makePlots
    plot_case_equilibrium_map(results.level, ...
        '6-DOF Level Equilibrium Map');
end

if saveCSV
    save_case_csv_equilibrium_map(results.level, ...
        'level_6dof_equilibrium_map_slow');
end

end

% =====================================================================
% LEVEL CASE
% =====================================================================

function out = run_level_case_6dof_equilibrium_map( ...
    db, Vlist, alt_list_m, ...
    PLA_min, PLA_max, DELTAZ_abs_bound_deg, ...
    showConsoleProgress, trim_tol, nConsecutiveNoEqStop)

nV = numel(Vlist);
nA = numel(alt_list_m);

PLA_req    = NaN(nA,nV);
DELTAZ_req = NaN(nA,nV);
alpha_req  = NaN(nA,nV);
theta_req  = NaN(nA,nV);
gamma_req  = NaN(nA,nV);
TW_req     = NaN(nA,nV);

equilibrium_ok = false(nA,nV);

failure_code   = zeros(nA,nV);
failure_reason = strings(nA,nV);
failure_reason(:) = "Not evaluated";
evaluated      = false(nA,nV);

exitflag_map   = NaN(nA,nV);
resnorm_map    = NaN(nA,nV);

x_trim_cell = cell(nA,nV);
u_trim_cell = cell(nA,nV);
trim_info   = cell(nA,nV);

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-8, ...
    'MaxFunctionEvaluations',1500, ...
    'MaxIterations',200);

for iv = 1:nV
    Vt = Vlist(iv);
    prev_sol = [];

    sawEquilibrium = false;
    consecutiveNoEq = 0;

    for ia = 1:nA
        alt_m = alt_list_m(ia);

        if isempty(prev_sol)
            y0 = default_trim_guess_level(Vt, alt_m, PLA_min, PLA_max);
        else
            y0 = prev_sol;
        end

        try
            [sol, info] = solve_trim_point_6dof_equilibrium( ...
                Vt, alt_m, 0.0, db, y0, opts, trim_tol, ...
                PLA_min, PLA_max, DELTAZ_abs_bound_deg);

            alpha_deg  = sol(1);
            theta_rad  = sol(2);
            PLA_deg    = sol(3);
            DELTAZ_deg = sol(4);

            [thisEq, failCode, failReason] = classify_trim_result_equilibrium_only( ...
                info, trim_tol, PLA_deg, PLA_min, PLA_max, DELTAZ_deg, DELTAZ_abs_bound_deg);

            alpha_req(ia,iv)  = alpha_deg;
            theta_req(ia,iv)  = rad2deg(theta_rad);
            gamma_req(ia,iv)  = info.gamma_deg;
            PLA_req(ia,iv)    = PLA_deg;
            DELTAZ_req(ia,iv) = DELTAZ_deg;
            TW_req(ia,iv)     = estimate_thrust_to_weight_ratio_from_trim(info, db);

            equilibrium_ok(ia,iv) = thisEq;
            evaluated(ia,iv)      = true;
            failure_code(ia,iv)   = failCode;
            failure_reason(ia,iv) = failReason;
            exitflag_map(ia,iv)   = info.exitflag;
            resnorm_map(ia,iv)    = info.resnorm;

            trim_info{ia,iv} = info;

            if thisEq
                x_trim_cell{ia,iv} = info.x_trim;
                u_trim_cell{ia,iv} = info.u_trim;
                prev_sol = sol;
                sawEquilibrium = true;
                consecutiveNoEq = 0;
            else
                if sawEquilibrium
                    consecutiveNoEq = consecutiveNoEq + 1;
                end
            end

            if showConsoleProgress
                fprintf(['LEVEL6-EQ | alt=%6.0f | Vt=%6.1f | eq=%d | ' ...
                         'alpha=%7.3f | theta=%7.3f | gamma=%7.3f | ' ...
                         'PLA=%9.3f | DELTAZ=%7.3f | TW=%8.4f | ' ...
                         'res=%9.3e | exit=%d | %s\n'], ...
                    alt_m, Vt, thisEq, ...
                    alpha_deg, rad2deg(theta_rad), info.gamma_deg, ...
                    PLA_deg, DELTAZ_deg, TW_req(ia,iv), ...
                    info.resnorm, info.exitflag, failReason);
            end

        catch ME
            evaluated(ia,iv)      = true;
            equilibrium_ok(ia,iv) = false;
            failure_code(ia,iv)   = 2;
            failure_reason(ia,iv) = "Solver/lookup error";
            exitflag_map(ia,iv)   = NaN;
            resnorm_map(ia,iv)    = NaN;

            trim_info{ia,iv} = struct('error_message', ME.message);

            if sawEquilibrium
                consecutiveNoEq = consecutiveNoEq + 1;
            end

            if showConsoleProgress
                fprintf('LEVEL6-EQ | alt=%6.0f | Vt=%6.1f | eq=0 | ERROR: %s\n', ...
                    alt_m, Vt, ME.message);
            end
        end

        if sawEquilibrium && consecutiveNoEq >= nConsecutiveNoEqStop
            if showConsoleProgress
                fprintf(['LEVEL6-EQ | alt=%6.0f | Vt=%6.1f | STOP after %d consecutive no-equilibrium points\n'], ...
                    alt_m, Vt, nConsecutiveNoEqStop);
            end
            break
        end
    end
end

out = struct();
out.caseName = 'level';
out.vy = 0;
out.Vlist = Vlist;
out.alt_list_m = alt_list_m;

out.PLA_req    = PLA_req;
out.DELTAZ_req = DELTAZ_req;
out.alpha_req  = alpha_req;
out.theta_req  = theta_req;
out.gamma_req  = gamma_req;
out.TW_req     = TW_req;

out.equilibrium_ok = equilibrium_ok;
out.failure_code   = failure_code;
out.failure_reason = failure_reason;
out.evaluated      = evaluated;
out.exitflag_map   = exitflag_map;
out.resnorm_map    = resnorm_map;

out.x_trim_cell = x_trim_cell;
out.u_trim_cell = u_trim_cell;
out.trim_info   = trim_info;

end

% =====================================================================
% CORE 6-DOF TRIM SOLVER
% =====================================================================

function [sol, info] = solve_trim_point_6dof_equilibrium( ...
    Vt, alt_m, gamma_target_rad, db, y0, opts, trim_tol, ...
    PLA_min, PLA_max, DELTAZ_abs_bound_deg)

obj = @(y) trim_residuals_6dof_safe( ...
    y, Vt, alt_m, gamma_target_rad, db, PLA_min, PLA_max, DELTAZ_abs_bound_deg);

[sol, fval, exitflag] = fsolve(obj, y0, opts);

[res, x_trim, u_trim, extra] = trim_residuals_6dof( ...
    sol, Vt, alt_m, gamma_target_rad, db, PLA_min, PLA_max, DELTAZ_abs_bound_deg);

resnorm = norm(res);

success = exitflag > 0 ...
    && all(isfinite(sol)) ...
    && all(isfinite(res)) ...
    && isfinite(resnorm) ...
    && (resnorm < trim_tol);

info = struct();
info.exitflag = exitflag;
info.fval = fval;
info.residual = res;
info.resnorm = resnorm;
info.success = success;
info.x_trim = x_trim;
info.u_trim = u_trim;
info.gamma_deg = rad2deg(extra.gamma_actual_rad);
info.gamma_err_deg = rad2deg(extra.gamma_actual_rad - gamma_target_rad);
info.dx = extra.dx;

end

function R = trim_residuals_6dof_safe( ...
    y, Vt, alt_m, gamma_target_rad, db, PLA_min, PLA_max, DELTAZ_abs_bound_deg)

try
    R = trim_residuals_6dof_only( ...
        y, Vt, alt_m, gamma_target_rad, db, PLA_min, PLA_max, DELTAZ_abs_bound_deg);

    if any(~isfinite(R))
        R = 1e6 * ones(4,1);
    end
catch
    R = 1e6 * ones(4,1);
end

end

function R = trim_residuals_6dof_only( ...
    y, Vt, alt_m, gamma_target_rad, db, PLA_min, PLA_max, DELTAZ_abs_bound_deg)

[R, ~, ~, ~] = trim_residuals_6dof( ...
    y, Vt, alt_m, gamma_target_rad, db, PLA_min, PLA_max, DELTAZ_abs_bound_deg);

end

function [R, x_trim, u_trim, extra] = trim_residuals_6dof( ...
    y, Vt, alt_m, gamma_target_rad, db, PLA_min, PLA_max, DELTAZ_abs_bound_deg)

% y = [alpha_deg; theta_rad; PLA_deg; DELTAZ_deg]

alpha_deg  = y(1);
theta_rad  = y(2);
PLA_deg    = y(3);
DELTAZ_deg = y(4);

if abs(alpha_deg) > 40 || ...
   abs(rad2deg(theta_rad)) > 50 || ...
   PLA_deg < PLA_min - 20 || PLA_deg > PLA_max || ...
   abs(DELTAZ_deg) > DELTAZ_abs_bound_deg

    R = 1e6 * ones(4,1);
    x_trim = nan(12,1);
    u_trim = struct('PLA',NaN,'DELTAZ',NaN,'DELTAW',NaN);
    extra = struct('gamma_actual_rad',NaN,'dx',nan(12,1));
    return
end

alpha_rad = deg2rad(alpha_deg);

u_b = Vt * cos(alpha_rad);
v_b = 0;
w_b = Vt * sin(alpha_rad);

x_trim = [ ...
    u_b; ...
    v_b; ...
    w_b; ...
    0; ...
    0; ...
    0; ...
    0; ...
    theta_rad; ...
    0; ...
    0; ...
    0; ...
   -alt_m];

u_trim = struct();
u_trim.PLA    = PLA_deg;
u_trim.DELTAZ = DELTAZ_deg;
u_trim.DELTAW = 0;

dx = ice_step_6dof(0, x_trim, u_trim, db);

gamma_actual_rad = gamma_from_state12_local(x_trim);

R = [ ...
    dx(1); ...
    dx(3); ...
    dx(5); ...
    gamma_actual_rad - gamma_target_rad ...
];

extra = struct();
extra.gamma_actual_rad = gamma_actual_rad;
extra.dx = dx;

end

% =====================================================================
% INITIAL GUESS
% =====================================================================

function y0 = default_trim_guess_level(Vt, alt_m, PLA_min, PLA_max)

alpha0_deg  = max(1, min(10, 2 + 0.015*(200 - Vt) + 0.0001*alt_m));
theta0_rad  = deg2rad(alpha0_deg);
PLA0_deg    = max(PLA_min + 5, min(PLA_max - 5, 50 + 0.004*alt_m + 0.04*(Vt - 100)));
DELTAZ0_deg = 0;

y0 = [alpha0_deg; theta0_rad; PLA0_deg; DELTAZ0_deg];

end

% =====================================================================
% SHARED HELPERS
% =====================================================================

function gamma_rad = gamma_from_state12_local(x12)

u = x12(1);
v = x12(2);
w = x12(3);

phi   = x12(7);
theta = x12(8);
psi   = x12(9);

cphi = cos(phi);   sphi = sin(phi);
cth  = cos(theta); sth  = sin(theta);
cpsi = cos(psi);   spsi = sin(psi);

R_b2n = [ cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
          cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
             -sth,                sphi*cth,                 cphi*cth      ];

Vn = R_b2n * [u; v; w];
gamma_rad = atan2(-Vn(3), sqrt(Vn(1)^2 + Vn(2)^2));

end

function [isEquilibrium, code, reason] = classify_trim_result_equilibrium_only( ...
    info, trim_tol, PLA_deg, PLA_min, PLA_max, DELTAZ_deg, DELTAZ_abs_bound_deg)

isEquilibrium = false;

if ~(isfinite(PLA_deg) && PLA_deg >= PLA_min - 20 && PLA_deg <= PLA_max)
    code = 6;
    reason = "PLA search bound";
    return
end

if ~(isfinite(DELTAZ_deg) && abs(DELTAZ_deg) <= DELTAZ_abs_bound_deg)
    code = 7;
    reason = "DELTAZ search bound";
    return
end

if info.exitflag <= 0
    code = 3;
    reason = "Non-converged exitflag";
    return
end

if ~(isfinite(info.resnorm) && info.resnorm < trim_tol)
    code = 4;
    reason = "Residual norm too large";
    return
end

if info.success
    isEquilibrium = true;
    code = 1;
    reason = "Equilibrium";
    return
end

code = 5;
reason = "Other failure";
end

function TW = estimate_thrust_to_weight_ratio_from_trim(info, db)

TW = NaN;

try
    x_trim = info.x_trim;
    u_trim = info.u_trim;

    if isempty(x_trim) || isempty(u_trim) || any(~isfinite(x_trim(1:12))) || ~isfield(u_trim,'PLA')
        return
    end

    u = x_trim(1);
    v = x_trim(2);
    w = x_trim(3);
    Vt = sqrt(u^2 + v^2 + w^2);

    eng = ice_block_engine(x_trim(12), Vt, u_trim.PLA, db);

    if ~isfield(eng,'FG') || ~isfinite(eng.FG)
        return
    end

    m = NaN;
    try, if isfield(db,'m'), m = db.m; end, catch, end
    try, if ~isfinite(m) && isfield(db,'mass'), m = db.mass; end, catch, end
    try, if ~isfinite(m) && isfield(db,'W'), m = db.W/9.81; end, catch, end

    if ~isfinite(m) || m <= 0
        return
    end

    TW = eng.FG / (m * 9.81);

catch
end

end

% =====================================================================
% CSV EXPORT
% =====================================================================

function save_case_csv_equilibrium_map(caseData, stem)

Vlist = caseData.Vlist;
alt_list_m = caseData.alt_list_m;
[Vgrid, AltGrid] = meshgrid(Vlist, alt_list_m);

T = table( ...
    Vgrid(:), AltGrid(:), ...
    caseData.equilibrium_ok(:), ...
    caseData.failure_code(:), ...
    cellstr(caseData.failure_reason(:)), ...
    caseData.exitflag_map(:), ...
    caseData.resnorm_map(:), ...
    caseData.alpha_req(:), ...
    caseData.theta_req(:), ...
    caseData.gamma_req(:), ...
    caseData.PLA_req(:), ...
    caseData.DELTAZ_req(:), ...
    caseData.TW_req(:), ...
    'VariableNames', { ...
    'Vt_mps','Altitude_m','EquilibriumOK','FailureCode','FailureReason', ...
    'Exitflag','Resnorm', ...
    'Alpha_deg','Theta_deg','Gamma_deg', ...
    'PLA_deg','DELTAZ_deg','TW_req'});

fname = matlab.lang.makeValidName(stem) + ".csv";
writetable(T, fname);
fprintf('Saved %s\n', fname);

end

% =====================================================================
% PLOTTING
% =====================================================================

function plot_case_equilibrium_map(caseData, figTitle)

Vlist = caseData.Vlist;
alt_list_m = caseData.alt_list_m;

figure('Color','w');
imagesc(Vlist, alt_list_m, caseData.equilibrium_ok);
set(gca,'YDir','normal')
colormap([1 1 1; 0.1 0.5 0.1])
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: Equilibrium Exists', figTitle))
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'No','Yes'};

figure('Color','w');
imagesc(Vlist, alt_list_m, caseData.PLA_req);
set(gca,'YDir','normal')
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: Required PLA', figTitle))
colorbar
colormap(gca, turbo)

figure('Color','w');
imagesc(Vlist, alt_list_m, caseData.DELTAZ_req);
set(gca,'YDir','normal')
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: Required DELTAZ', figTitle))
colorbar
colormap(gca, turbo)

figure('Color','w');
imagesc(Vlist, alt_list_m, caseData.TW_req);
set(gca,'YDir','normal')
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: Approximate thrust-to-weight ratio', figTitle))
colorbar
colormap(gca, turbo)

end