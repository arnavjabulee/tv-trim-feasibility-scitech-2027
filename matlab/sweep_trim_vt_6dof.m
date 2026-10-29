function results = sweep_trim_vt_6dof(mode)
% SWEEP_TRIM_VT_6DOF
%
% Finds actual 6-DOF trim points using ice_step_6dof + fsolve.
%
% Modes:
%   'all'   -> run steady level + steady climb/descent
%   'level' -> run steady level only
%   'climb' -> run steady climb/descent only
%
% This version sweeps multiple DELTAZ limits and labels outputs.
%
% 12-state order:
%   x = [u v w p q r phi theta psi xN yE zD]

if nargin < 1
    mode = 'level';
end

mode = lower(string(mode));

db = ice_init();

% ---------------- USER OPTIONS ----------------
Vlist = 76:2:560;              % m/s
alt_list_m = 0:100:20000;      % m
vy_list = -20:2:20;            % m/s

DELTAZ_lim_list = [15 20 25];  % deg
PLA_min = min(db.eng.PLA);
PLA_max = max(db.eng.PLA);

showConsoleProgress = true;

% For debugging, set to inf
nConsecutiveInfeasibleStop = 8;

trim_tol = 1e-5;
% ---------------------------------------------

results = struct();
results.deltaz_limits_deg = DELTAZ_lim_list;

for iLim = 1:numel(DELTAZ_lim_list)
    DELTAZ_lim = DELTAZ_lim_list(iLim);
    limField = matlab.lang.makeValidName(sprintf('deltaz_%ddeg', DELTAZ_lim));

    fprintf('\n============================================================\n');
    fprintf('Running TVC trim sweep for DELTAZ limit = %g deg\n', DELTAZ_lim);
    fprintf('============================================================\n');

    results.(limField) = struct();
    results.(limField).DELTAZ_lim_deg = DELTAZ_lim;

    if mode == "all" || mode == "level"
        results.(limField).level = run_level_case_6dof( ...
            db, Vlist, alt_list_m, PLA_min, PLA_max, DELTAZ_lim, ...
            showConsoleProgress, nConsecutiveInfeasibleStop, trim_tol);

        plot_case(results.(limField).level, db, ...
            sprintf('6-DOF Level Trim (DELTAZ limit = %g deg)', DELTAZ_lim), ...
            DELTAZ_lim);

        save_case_csv(results.(limField).level, ...
            sprintf('level_6dof_deltaz_%gdeg', DELTAZ_lim));
    end

    if mode == "all" || mode == "climb"
        climbCases = cell(numel(vy_list),1);

        for k = 1:numel(vy_list)
            vy = vy_list(k);

            climbCases{k} = run_climb_case_6dof( ...
                db, Vlist, alt_list_m, vy, PLA_min, PLA_max, DELTAZ_lim, ...
                showConsoleProgress, nConsecutiveInfeasibleStop, trim_tol);

            plot_case(climbCases{k}, db, ...
                sprintf('6-DOF Steady Climb/Descent Trim (v_y = %.1f m/s, DELTAZ limit = %g deg)', ...
                vy, DELTAZ_lim), ...
                DELTAZ_lim);

            save_case_csv(climbCases{k}, ...
                sprintf('climb6dof_vy_%+.1f_deltaz_%gdeg', vy, DELTAZ_lim));
        end

        results.(limField).climb = climbCases;
    end
end

end

% =====================================================================
% LEVEL CASE
% =====================================================================

function out = run_level_case_6dof( ...
    db, Vlist, alt_list_m, PLA_min, PLA_max, DELTAZ_lim, ...
    showConsoleProgress, nConsecutiveInfeasibleStop, trim_tol)

nV = numel(Vlist);
nA = numel(alt_list_m);

PLA_req    = NaN(nA,nV);
DELTAZ_req = NaN(nA,nV);
alpha_req  = NaN(nA,nV);
theta_req  = NaN(nA,nV);
gamma_req  = NaN(nA,nV);
T_req      = NaN(nA,nV);
feasible   = false(nA,nV);

failure_code   = zeros(nA,nV);
failure_reason = strings(nA,nV);
failure_reason(:) = "Not evaluated";
evaluated      = false(nA,nV);

x_trim_cell = cell(nA,nV);
u_trim_cell = cell(nA,nV);
trim_info   = cell(nA,nV);

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-8, ...
    'MaxFunctionEvaluations',1200, ...
    'MaxIterations',150);

for iv = 1:nV
    Vt = Vlist(iv);

    sawFeasible = false;
    consecutiveInfeasible = 0;
    prev_sol = [];

    for ia = 1:nA
        alt_m = alt_list_m(ia);

        if isempty(prev_sol)
            y0 = default_trim_guess_level(Vt, alt_m, PLA_min, PLA_max);
        else
            y0 = prev_sol;
        end

        try
            [sol, info] = solve_trim_point_6dof( ...
                Vt, alt_m, 0.0, db, y0, opts, trim_tol);

            alpha_deg  = sol(1);
            theta_rad  = sol(2);
            PLA_deg    = sol(3);
            DELTAZ_deg = sol(4);

            [thisFeasible, failCode, failReason] = classify_trim_result_tvc( ...
                info, PLA_deg, DELTAZ_deg, PLA_min, PLA_max, DELTAZ_lim, trim_tol);

            alpha_req(ia,iv)  = alpha_deg;
            theta_req(ia,iv)  = rad2deg(theta_rad);
            gamma_req(ia,iv)  = info.gamma_deg;
            PLA_req(ia,iv)    = PLA_deg;
            DELTAZ_req(ia,iv) = DELTAZ_deg;
            T_req(ia,iv)      = NaN;
            feasible(ia,iv)   = thisFeasible;

            evaluated(ia,iv)      = true;
            failure_code(ia,iv)   = failCode;
            failure_reason(ia,iv) = failReason;

            trim_info{ia,iv} = info;

            if thisFeasible
                x_trim_cell{ia,iv} = info.x_trim;
                u_trim_cell{ia,iv} = info.u_trim;
                prev_sol = sol;
            end

            if showConsoleProgress
                fprintf(['LEVEL6 | DELTAZlim=%5.1f | alt=%6.0f | Vt=%6.1f | feas=%d | ' ...
                         'alpha=%7.3f | theta=%7.3f | gamma=%7.3f | ' ...
                         'PLA=%7.3f | DELTAZ=%7.3f | res=%9.3e | exit=%d\n'], ...
                    DELTAZ_lim, alt_m, Vt, thisFeasible, ...
                    alpha_deg, rad2deg(theta_rad), info.gamma_deg, ...
                    PLA_deg, DELTAZ_deg, info.resnorm, info.exitflag);

                fprintf('         residuals: du=% .3e | dw=% .3e | dq=% .3e | dgamma=% .3e\n', ...
                    info.residual(1), info.residual(2), info.residual(3), info.residual(4));

                if ~thisFeasible
                    fprintf('         fail reason: %s\n', failReason);
                end
            end

        catch ME
            thisFeasible = false;

            evaluated(ia,iv)      = true;
            feasible(ia,iv)       = false;
            failure_code(ia,iv)   = 2;
            failure_reason(ia,iv) = "Solver/lookup error";

            trim_info{ia,iv} = struct('error_message', ME.message);

            if showConsoleProgress
                fprintf('LEVEL6 | DELTAZlim=%5.1f | alt=%6.0f | Vt=%6.1f | feas=0 | ERROR: %s\n', ...
                    DELTAZ_lim, alt_m, Vt, ME.message);
            end
        end

        if thisFeasible
            sawFeasible = true;
            consecutiveInfeasible = 0;
        else
            if sawFeasible
                consecutiveInfeasible = consecutiveInfeasible + 1;
                if consecutiveInfeasible >= nConsecutiveInfeasibleStop
                    break
                end
            end
        end
    end
end

out = struct();
out.caseName = 'level';
out.vy = 0;
out.DELTAZ_lim_deg = DELTAZ_lim;
out.Vlist = Vlist;
out.alt_list_m = alt_list_m;
out.PLA_req = PLA_req;
out.DELTAZ_req = DELTAZ_req;
out.alpha_req = alpha_req;
out.theta_req = theta_req;
out.gamma_req = gamma_req;
out.T_req = T_req;
out.feasible = feasible;
out.failure_code = failure_code;
out.failure_reason = failure_reason;
out.evaluated = evaluated;
out.x_trim_cell = x_trim_cell;
out.u_trim_cell = u_trim_cell;
out.trim_info = trim_info;
out.border_fail_mask = find_border_failures(feasible, evaluated);
report_border_failures(out);

end

% =====================================================================
% CLIMB / DESCENT CASE
% =====================================================================

function out = run_climb_case_6dof( ...
    db, Vlist, alt_list_m, vy, PLA_min, PLA_max, DELTAZ_lim, ...
    showConsoleProgress, nConsecutiveInfeasibleStop, trim_tol)

nV = numel(Vlist);
nA = numel(alt_list_m);

PLA_req    = NaN(nA,nV);
DELTAZ_req = NaN(nA,nV);
alpha_req  = NaN(nA,nV);
theta_req  = NaN(nA,nV);
gamma_req  = NaN(nA,nV);
T_req      = NaN(nA,nV);
feasible   = false(nA,nV);

failure_code   = zeros(nA,nV);
failure_reason = strings(nA,nV);
failure_reason(:) = "Not evaluated";
evaluated      = false(nA,nV);

x_trim_cell = cell(nA,nV);
u_trim_cell = cell(nA,nV);
trim_info   = cell(nA,nV);

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-8, ...
    'MaxFunctionEvaluations',1200, ...
    'MaxIterations',150);

for iv = 1:nV
    Vt = Vlist(iv);

    sawFeasible = false;
    consecutiveInfeasible = 0;
    prev_sol = [];

    if abs(vy) >= Vt
        continue
    end

    gamma_target_rad = asin(vy / Vt);

    for ia = 1:nA
        alt_m = alt_list_m(ia);

        if isempty(prev_sol)
            y0 = default_trim_guess_climb(Vt, alt_m, gamma_target_rad, PLA_min, PLA_max);
        else
            y0 = prev_sol;
        end

        try
            [sol, info] = solve_trim_point_6dof( ...
                Vt, alt_m, gamma_target_rad, db, y0, opts, trim_tol);

            alpha_deg  = sol(1);
            theta_rad  = sol(2);
            PLA_deg    = sol(3);
            DELTAZ_deg = sol(4);

            [thisFeasible, failCode, failReason] = classify_trim_result_tvc( ...
                info, PLA_deg, DELTAZ_deg, PLA_min, PLA_max, DELTAZ_lim, trim_tol);

            alpha_req(ia,iv)  = alpha_deg;
            theta_req(ia,iv)  = rad2deg(theta_rad);
            gamma_req(ia,iv)  = info.gamma_deg;
            PLA_req(ia,iv)    = PLA_deg;
            DELTAZ_req(ia,iv) = DELTAZ_deg;
            T_req(ia,iv)      = NaN;
            feasible(ia,iv)   = thisFeasible;

            evaluated(ia,iv)      = true;
            failure_code(ia,iv)   = failCode;
            failure_reason(ia,iv) = failReason;

            trim_info{ia,iv} = info;

            if thisFeasible
                x_trim_cell{ia,iv} = info.x_trim;
                u_trim_cell{ia,iv} = info.u_trim;
                prev_sol = sol;
            end

            if showConsoleProgress
                fprintf(['CLIMB6 | DELTAZlim=%5.1f | vy=%6.1f | alt=%6.0f | Vt=%6.1f | feas=%d | ' ...
                         'alpha=%7.3f | theta=%7.3f | gamma=%7.3f | ' ...
                         'PLA=%7.3f | DELTAZ=%7.3f | res=%9.3e | exit=%d\n'], ...
                    DELTAZ_lim, vy, alt_m, Vt, thisFeasible, ...
                    alpha_deg, rad2deg(theta_rad), info.gamma_deg, ...
                    PLA_deg, DELTAZ_deg, info.resnorm, info.exitflag);

                fprintf('         residuals: du=% .3e | dw=% .3e | dq=% .3e | dgamma=% .3e\n', ...
                    info.residual(1), info.residual(2), info.residual(3), info.residual(4));

                if ~thisFeasible
                    fprintf('         fail reason: %s\n', failReason);
                end
            end

        catch ME
            thisFeasible = false;

            evaluated(ia,iv)      = true;
            feasible(ia,iv)       = false;
            failure_code(ia,iv)   = 2;
            failure_reason(ia,iv) = "Solver/lookup error";

            trim_info{ia,iv} = struct('error_message', ME.message);

            if showConsoleProgress
                fprintf('CLIMB6 | DELTAZlim=%5.1f | vy=%6.1f | alt=%6.0f | Vt=%6.1f | feas=0 | ERROR: %s\n', ...
                    DELTAZ_lim, vy, alt_m, Vt, ME.message);
            end
        end

        if thisFeasible
            sawFeasible = true;
            consecutiveInfeasible = 0;
        else
            if sawFeasible
                consecutiveInfeasible = consecutiveInfeasible + 1;
                if consecutiveInfeasible >= nConsecutiveInfeasibleStop
                    break
                end
            end
        end
    end
end

out = struct();
out.caseName = 'climb';
out.vy = vy;
out.DELTAZ_lim_deg = DELTAZ_lim;
out.Vlist = Vlist;
out.alt_list_m = alt_list_m;
out.PLA_req = PLA_req;
out.DELTAZ_req = DELTAZ_req;
out.alpha_req = alpha_req;
out.theta_req = theta_req;
out.gamma_req = gamma_req;
out.T_req = T_req;
out.feasible = feasible;
out.failure_code = failure_code;
out.failure_reason = failure_reason;
out.evaluated = evaluated;
out.x_trim_cell = x_trim_cell;
out.u_trim_cell = u_trim_cell;
out.trim_info = trim_info;
out.border_fail_mask = find_border_failures(feasible, evaluated);
report_border_failures(out);

end

% =====================================================================
% CORE 6-DOF TRIM SOLVER
% =====================================================================

function [sol, info] = solve_trim_point_6dof(Vt, alt_m, gamma_target_rad, db, y0, opts, trim_tol)

obj = @(y) trim_residuals_6dof_safe(y, Vt, alt_m, gamma_target_rad, db);

[sol, fval, exitflag] = fsolve(obj, y0, opts);

[res, x_trim, u_trim, extra] = trim_residuals_6dof(sol, Vt, alt_m, gamma_target_rad, db);
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

function R = trim_residuals_6dof_safe(y, Vt, alt_m, gamma_target_rad, db)

try
    R = trim_residuals_6dof_only(y, Vt, alt_m, gamma_target_rad, db);

    if any(~isfinite(R))
        R = 1e6 * ones(4,1);
    end
catch
    R = 1e6 * ones(4,1);
end

end

function R = trim_residuals_6dof_only(y, Vt, alt_m, gamma_target_rad, db)
[R, ~, ~, ~] = trim_residuals_6dof(y, Vt, alt_m, gamma_target_rad, db);
end

function [R, x_trim, u_trim, extra] = trim_residuals_6dof(y, Vt, alt_m, gamma_target_rad, db)

% Unknowns:
% y = [alpha_deg; theta_rad; PLA_deg; DELTAZ_deg]

alpha_deg  = y(1);
theta_rad  = y(2);
PLA_deg    = y(3);
DELTAZ_deg = y(4);

if abs(alpha_deg) > 40 || ...
   abs(rad2deg(theta_rad)) > 50 || ...
   PLA_deg < -10 || PLA_deg > 150 || ...
   abs(DELTAZ_deg) > 50

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
% INITIAL GUESSES
% =====================================================================

function y0 = default_trim_guess_level(Vt, alt_m, PLA_min, PLA_max)

alpha0_deg = max(1, min(10, 2 + 0.015*(200 - Vt) + 0.0001*alt_m));
theta0_rad = deg2rad(alpha0_deg);
PLA0_deg   = max(PLA_min + 5, min(PLA_max - 5, 50 + 0.004*alt_m + 0.04*(Vt - 100)));
DELTAZ0_deg = 0;

y0 = [alpha0_deg; theta0_rad; PLA0_deg; DELTAZ0_deg];

end

function y0 = default_trim_guess_climb(Vt, alt_m, gamma_target_rad, PLA_min, PLA_max)

gamma_deg = rad2deg(gamma_target_rad);
alpha0_deg = max(1, min(12, 3 + 0.015*(200 - Vt) + 0.0001*alt_m));
theta0_rad = deg2rad(alpha0_deg + gamma_deg);
PLA0_deg   = max(PLA_min + 5, min(PLA_max - 5, 55 + 0.004*alt_m + 0.04*(Vt - 100) + 0.8*abs(gamma_deg)));
DELTAZ0_deg = 0;

y0 = [alpha0_deg; theta0_rad; PLA0_deg; DELTAZ0_deg];

end

% =====================================================================
% GAMMA FROM 12-STATE
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

% =====================================================================
% FAILURE ANALYSIS HELPERS
% =====================================================================

function [isFeasible, code, reason] = classify_trim_result_tvc( ...
    info, PLA_deg, DELTAZ_deg, PLA_min, PLA_max, DELTAZ_lim, trim_tol)

% 1 feasible
% 2 solver/lookup error
% 3 exitflag/non-convergence
% 4 PLA limit
% 5 DELTAZ limit
% 6 residual norm tolerance
% 7 other

isFeasible = false;

if ~(isfinite(PLA_deg) && PLA_deg >= PLA_min && PLA_deg <= PLA_max)
    code = 4;
    reason = "PLA limit";
    return
end

if ~(isfinite(DELTAZ_deg) && abs(DELTAZ_deg) <= DELTAZ_lim)
    code = 5;
    reason = "DELTAZ limit";
    return
end

if info.exitflag <= 0
    code = 3;
    reason = "Non-converged exitflag";
    return
end

if ~(isfinite(info.resnorm) && info.resnorm < trim_tol)
    code = 6;
    reason = "Residual norm tolerance";
    return
end

if info.success
    isFeasible = true;
    code = 1;
    reason = "Feasible";
    return
end

code = 7;
reason = "Other failure";
end

function border_fail_mask = find_border_failures(feasible, evaluated)
[nA, nV] = size(feasible);
border_fail_mask = false(nA,nV);

for ia = 1:nA
    for iv = 1:nV
        if evaluated(ia,iv) && ~feasible(ia,iv)
            r1 = max(1, ia-1); r2 = min(nA, ia+1);
            c1 = max(1, iv-1); c2 = min(nV, iv+1);
            nbr = feasible(r1:r2, c1:c2);
            if any(nbr(:))
                border_fail_mask(ia,iv) = true;
            end
        end
    end
end
end

function report_border_failures(caseData)
mask = caseData.border_fail_mask;
if ~any(mask(:))
    fprintf('\nNo border failures found for case "%s" (DELTAZ limit = %.1f deg).\n', ...
        caseData.caseName, caseData.DELTAZ_lim_deg);
    return
end

codes = caseData.failure_code(mask);
reasons = caseData.failure_reason(mask);

fprintf('\nBorder failure summary for case "%s" (DELTAZ limit = %.1f deg):\n', ...
    caseData.caseName, caseData.DELTAZ_lim_deg);

uc = unique(codes(:))';
for k = uc
    if k <= 0, continue; end
    n = sum(codes == k);
    idx = find(codes == k, 1, 'first');
    fprintf('  code %d : %4d cells : %s\n', k, n, reasons(idx));
end

[row, col] = find(mask);
fprintf('  first few border failures:\n');
nshow = min(15, numel(row));
for i = 1:nshow
    ia = row(i);
    iv = col(i);
    fprintf('    Alt=%7.1f m, Vt=%7.1f m/s, reason=%s\n', ...
        caseData.alt_list_m(ia), caseData.Vlist(iv), caseData.failure_reason(ia,iv));
end
fprintf('\n');
end

% =====================================================================
% CSV EXPORT
% =====================================================================

function save_case_csv(caseData, stem)

Vlist = caseData.Vlist;
alt_list_m = caseData.alt_list_m;
[Vgrid, AltGrid] = meshgrid(Vlist, alt_list_m);

T = table( ...
    Vgrid(:), AltGrid(:), ...
    caseData.feasible(:), ...
    caseData.failure_code(:), ...
    cellstr(caseData.failure_reason(:)), ...
    caseData.alpha_req(:), ...
    caseData.theta_req(:), ...
    caseData.gamma_req(:), ...
    caseData.PLA_req(:), ...
    caseData.DELTAZ_req(:), ...
    caseData.T_req(:), ...
    'VariableNames', { ...
    'Vt_mps','Altitude_m','Feasible','FailureCode','FailureReason', ...
    'Alpha_deg','Theta_deg','Gamma_deg', ...
    'PLA_deg','DELTAZ_deg','Treq_N'});

fname = matlab.lang.makeValidName(stem) + ".csv";
writetable(T, fname);
fprintf('Saved %s\n', fname);

end

% =====================================================================
% PLOTTING
% =====================================================================

function plot_case(caseData, db, figTitle, DELTAZ_lim)

Vlist = caseData.Vlist;
alt_list_m = caseData.alt_list_m;
PLA_req = caseData.PLA_req;
DELTAZ_req = caseData.DELTAZ_req;
feasible = caseData.feasible;
failure_code = caseData.failure_code;
border_fail_mask = caseData.border_fail_mask;

[Vgrid, AltGrid] = meshgrid(Vlist, alt_list_m);

PLA_min = min(db.eng.PLA);
PLA_max = max(db.eng.PLA);

% ---------- PLA ----------
figure('Color','w');
hold on

PLA_feas = PLA_req;
PLA_infeas = PLA_req;
PLA_feas(~feasible) = NaN;
PLA_infeas(feasible) = NaN;

s1 = surf(Vgrid, AltGrid, PLA_feas, PLA_feas, ...
    'EdgeColor','none', ...
    'FaceColor','interp');

colormap(gca, turbo)
cb1 = colorbar;
cb1.Label.String = 'Required PLA (deg)';

PLA_violation = zeros(size(PLA_req));
PLA_violation(PLA_req > PLA_max) = PLA_req(PLA_req > PLA_max) - PLA_max;
PLA_violation(PLA_req < PLA_min) = PLA_min - PLA_req(PLA_req < PLA_min);
PLA_violation(feasible) = NaN;

redmap = zeros([size(PLA_req), 3]);
if any(~isnan(PLA_violation(:)))
    vmax = max(PLA_violation(:), [], 'omitnan');
    if vmax == 0, vmax = 1; end
    vnorm = PLA_violation / vmax;
    redmap(:,:,1) = 1.0;
    redmap(:,:,2) = 0.9 - 0.7*vnorm;
    redmap(:,:,3) = 0.9 - 0.7*vnorm;
else
    redmap(:,:,1) = 1.0;
    redmap(:,:,2) = 0.9;
    redmap(:,:,3) = 0.9;
end

s2 = surf(Vgrid, AltGrid, PLA_infeas, redmap, ...
    'EdgeColor','none', ...
    'FaceColor','interp');

s3 = surf(Vgrid, AltGrid, PLA_max*ones(size(PLA_req)), ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'FaceColor',[0 0 0]);

s4 = surf(Vgrid, AltGrid, PLA_min*ones(size(PLA_req)), ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'FaceColor',[0 0 0]);

grid on
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
zlabel('PLA (deg)')
title(sprintf('%s: Required PLA', figTitle))
view(45,30)

legend([s1 s2 s3 s4], ...
    {'Feasible PLA', ...
     'Out-of-bounds PLA', ...
     'PLA upper limit', ...
     'PLA lower limit'}, ...
    'Location','best');

% ---------- DELTAZ ----------
figure('Color','w');
hold on

DZ_feas = DELTAZ_req;
DZ_infeas = DELTAZ_req;
DZ_feas(~feasible) = NaN;
DZ_infeas(feasible) = NaN;

s1 = surf(Vgrid, AltGrid, DZ_feas, DZ_feas, ...
    'EdgeColor','none', ...
    'FaceColor','interp');

colormap(gca, turbo)
cb2 = colorbar;
cb2.Label.String = 'Required \Delta Z (deg)';

DZ_violation = zeros(size(DELTAZ_req));
DZ_violation(DELTAZ_req > DELTAZ_lim)  = DELTAZ_req(DELTAZ_req > DELTAZ_lim) - DELTAZ_lim;
DZ_violation(DELTAZ_req < -DELTAZ_lim) = -DELTAZ_lim - DELTAZ_req(DELTAZ_req < -DELTAZ_lim);
DZ_violation(feasible) = NaN;

redmap = zeros([size(DELTAZ_req), 3]);
if any(~isnan(DZ_violation(:)))
    vmax = max(DZ_violation(:), [], 'omitnan');
    if vmax == 0, vmax = 1; end
    vnorm = DZ_violation / vmax;
    redmap(:,:,1) = 1.0;
    redmap(:,:,2) = 0.9 - 0.7*vnorm;
    redmap(:,:,3) = 0.9 - 0.7*vnorm;
else
    redmap(:,:,1) = 1.0;
    redmap(:,:,2) = 0.9;
    redmap(:,:,3) = 0.9;
end

s2 = surf(Vgrid, AltGrid, DZ_infeas, redmap, ...
    'EdgeColor','none', ...
    'FaceColor','interp');

s3 = surf(Vgrid, AltGrid, DELTAZ_lim*ones(size(DZ_feas)), ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'FaceColor',[0 0 0]);

s4 = surf(Vgrid, AltGrid, -DELTAZ_lim*ones(size(DZ_feas)), ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'FaceColor',[0 0 0]);

grid on
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
zlabel('\Delta Z (deg)')
title(sprintf('%s: Required DELTAZ', figTitle))
view(45,30)

legend([s1 s2 s3 s4], ...
    {'Feasible \Delta Z', ...
     'Out-of-bounds \Delta Z', ...
     '\Delta Z upper limit', ...
     '\Delta Z lower limit'}, ...
    'Location','best');

% ---------- 2D envelope ----------
figure('Color','w');
imagesc(Vlist, alt_list_m, feasible);
set(gca,'YDir','normal')
colormap([1 1 1; 0 0.6 0])
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: Feasibility Envelope', figTitle))
cb = colorbar;
cb.Ticks = [0 1];
cb.TickLabels = {'Infeasible','Feasible'};

% ---------- Border failure reasons ----------
display_code = zeros(size(failure_code));
display_code(feasible) = 1;
display_code(border_fail_mask) = failure_code(border_fail_mask);

figure('Color','w');
imagesc(Vlist, alt_list_m, display_code);
set(gca,'YDir','normal')

% 0 white
% 1 green feasible
% 2 black solver
% 3 grey exitflag
% 4 red PLA
% 5 orange DELTAZ
% 6 blue residual norm
% 7 brown other
cmap = [ ...
    1.00 1.00 1.00;   % 0 not shown / non-border infeasible
    0.00 0.60 0.00;   % 1 feasible
    0.00 0.00 0.00;   % 2 solver
    0.45 0.45 0.45;   % 3 exitflag
    0.85 0.10 0.10;   % 4 PLA
    1.00 0.50 0.00;   % 5 DELTAZ
    0.10 0.30 0.90;   % 6 residual norm
    0.55 0.35 0.15];  % 7 other

colormap(gca, cmap)
caxis([0 7])
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: Border Failure Reasons', figTitle))

cb = colorbar;
cb.Ticks = 0:7;
cb.TickLabels = { ...
    'Not shown', ...
    'Feasible', ...
    'Solver error', ...
    'Exitflag', ...
    'PLA limit', ...
    'DELTAZ limit', ...
    'Residual norm', ...
    'Other'};
end