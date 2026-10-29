function results = sweep_trim_vt_6dof_aero(mode)
% SWEEP_TRIM_VT_6DOF_AERO
%
% Aerodynamic-controls trim sweep using:
%   y = [alpha_deg; theta_rad; PLA_deg; PFLAP_deg; SPLRL_deg]

if nargin < 1
    mode = 'level';
end

mode = lower(string(mode));
db = ice_init_aeroctrl();

% ---------------- USER OPTIONS ----------------
Vlist = 80:2:550;               % m/s
alt_list_m = 0:100:20000;       % m
vy_list = -20:2:20;             % m/s

PFLAP_lim = max(abs(db.aero.ctrl.PFLAP));
SPLRL_lim = max(abs(db.aero.ctrl.SPLRL2));
PLA_min = min(db.eng.PLA);
PLA_max = max(db.eng.PLA);

showConsoleProgress = true;
nConsecutiveInfeasibleStop = 8;   % set inf if you want full sweep always

reg = struct();
reg.w_pflap = 1e-3;
reg.w_splrl = 1e-3;
% ---------------------------------------------

results = struct();

if mode == "all" || mode == "level"
    results.level = run_level_case_6dof_aero( ...
        db, Vlist, alt_list_m, PLA_min, PLA_max, PFLAP_lim, SPLRL_lim, ...
        showConsoleProgress, nConsecutiveInfeasibleStop, reg);

    plot_case_aero(results.level, db, '6-DOF Level Trim (Aero Ctrl)');
    save_case_csv_aero(results.level, 'level_6dof_aero');
end

if mode == "all" || mode == "climb"
    climbCases = cell(numel(vy_list),1);

    for k = 1:numel(vy_list)
        vy = vy_list(k);

        climbCases{k} = run_climb_case_6dof_aero( ...
            db, Vlist, alt_list_m, vy, PLA_min, PLA_max, PFLAP_lim, SPLRL_lim, ...
            showConsoleProgress, nConsecutiveInfeasibleStop, reg);

        plot_case_aero(climbCases{k}, db, ...
            sprintf('6-DOF Steady Climb/Descent Trim (Aero Ctrl, v_y = %.1f m/s)', vy));

        save_case_csv_aero(climbCases{k}, sprintf('climb6dof_aero_vy_%+.1f', vy));
    end

    results.climb = climbCases;
end
end

% =====================================================================
% LEVEL CASE
% =====================================================================

function out = run_level_case_6dof_aero( ...
    db, Vlist, alt_list_m, PLA_min, PLA_max, PFLAP_lim, SPLRL_lim, ...
    showConsoleProgress, nConsecutiveInfeasibleStop, reg)

nV = numel(Vlist);
nA = numel(alt_list_m);

PLA_req   = NaN(nA,nV);
PFLAP_req = NaN(nA,nV);
SPLRL_req = NaN(nA,nV);
alpha_req = NaN(nA,nV);
theta_req = NaN(nA,nV);
gamma_req = NaN(nA,nV);
feasible  = false(nA,nV);

failure_code   = zeros(nA,nV);              % 0 = not evaluated
failure_reason = strings(nA,nV);
failure_reason(:) = "Not evaluated";
evaluated      = false(nA,nV);

x_trim_cell = cell(nA,nV);
u_trim_cell = cell(nA,nV);
trim_info   = cell(nA,nV);

opts = optimoptions('lsqnonlin', ...
    'Display','off', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10, ...
    'MaxFunctionEvaluations',2500, ...
    'MaxIterations',300);

lb = [-40; deg2rad(-50); PLA_min; -PFLAP_lim; -SPLRL_lim];
ub = [ 40; deg2rad( 50); PLA_max;  PFLAP_lim;  SPLRL_lim];

for iv = 1:nV
    Vt = Vlist(iv);

    sawFeasible = false;
    consecutiveInfeasible = 0;
    prev_sol = [];

    for ia = 1:nA
        alt_m = alt_list_m(ia);

        if isempty(prev_sol)
            y0 = default_trim_guess_level_aero(Vt, alt_m, PLA_min, PLA_max);
        else
            y0 = prev_sol;
        end
        y0 = min(max(y0, lb), ub);

        try
            [sol, info] = solve_trim_point_6dof_aero( ...
                Vt, alt_m, 0.0, db, y0, lb, ub, opts, reg);

            alpha_deg = sol(1);
            theta_rad = sol(2);
            PLA_deg   = sol(3);
            PFLAP_deg = sol(4);
            SPLRL_deg = sol(5);

            [thisFeasible, failCode, failReason] = classify_trim_result( ...
                info, PLA_deg, PFLAP_deg, SPLRL_deg, ...
                PLA_min, PLA_max, PFLAP_lim, SPLRL_lim);

            alpha_req(ia,iv) = alpha_deg;
            theta_req(ia,iv) = rad2deg(theta_rad);
            gamma_req(ia,iv) = info.gamma_deg;
            PLA_req(ia,iv)   = PLA_deg;
            PFLAP_req(ia,iv) = PFLAP_deg;
            SPLRL_req(ia,iv) = SPLRL_deg;
            feasible(ia,iv)  = thisFeasible;

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
                fprintf(['AERO-LVL | alt=%6.0f | Vt=%6.1f | feas=%d | ' ...
                         'alpha=%7.3f | theta=%7.3f | gamma=%7.3f | ' ...
                         'PLA=%7.3f | PFLAP=%7.3f | SPLRL=%7.3f | res=%9.3e | exit=%d\n'], ...
                    alt_m, Vt, thisFeasible, alpha_deg, rad2deg(theta_rad), ...
                    info.gamma_deg, PLA_deg, PFLAP_deg, SPLRL_deg, ...
                    info.resnorm_phys, info.exitflag);

                fprintf('           residuals: du=% .3e | dw=% .3e | dq=% .3e | dgamma=% .3e\n', ...
                    info.residual_phys(1), info.residual_phys(2), ...
                    info.residual_phys(3), info.residual_phys(4));

                if ~thisFeasible
                    fprintf('           fail reason: %s\n', failReason);
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
                fprintf('AERO-LVL | alt=%6.0f | Vt=%6.1f | feas=0 | ERROR: %s\n', ...
                    alt_m, Vt, ME.message);
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
out.Vlist = Vlist;
out.alt_list_m = alt_list_m;
out.PLA_req = PLA_req;
out.PFLAP_req = PFLAP_req;
out.SPLRL_req = SPLRL_req;
out.alpha_req = alpha_req;
out.theta_req = theta_req;
out.gamma_req = gamma_req;
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

function out = run_climb_case_6dof_aero( ...
    db, Vlist, alt_list_m, vy, PLA_min, PLA_max, PFLAP_lim, SPLRL_lim, ...
    showConsoleProgress, nConsecutiveInfeasibleStop, reg)

nV = numel(Vlist);
nA = numel(alt_list_m);

PLA_req   = NaN(nA,nV);
PFLAP_req = NaN(nA,nV);
SPLRL_req = NaN(nA,nV);
alpha_req = NaN(nA,nV);
theta_req = NaN(nA,nV);
gamma_req = NaN(nA,nV);
feasible  = false(nA,nV);

failure_code   = zeros(nA,nV);
failure_reason = strings(nA,nV);
failure_reason(:) = "Not evaluated";
evaluated      = false(nA,nV);

x_trim_cell = cell(nA,nV);
u_trim_cell = cell(nA,nV);
trim_info   = cell(nA,nV);

opts = optimoptions('lsqnonlin', ...
    'Display','off', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10, ...
    'MaxFunctionEvaluations',2500, ...
    'MaxIterations',300);

lb = [-40; deg2rad(-50); PLA_min; -PFLAP_lim; -SPLRL_lim];
ub = [ 40; deg2rad( 50); PLA_max;  PFLAP_lim;  SPLRL_lim];

for iv = 1:nV
    Vt = Vlist(iv);
    if abs(vy) >= Vt
        continue
    end

    gamma_target_rad = asin(vy / Vt);

    sawFeasible = false;
    consecutiveInfeasible = 0;
    prev_sol = [];

    for ia = 1:nA
        alt_m = alt_list_m(ia);

        if isempty(prev_sol)
            y0 = default_trim_guess_climb_aero(Vt, alt_m, gamma_target_rad, PLA_min, PLA_max);
        else
            y0 = prev_sol;
        end
        y0 = min(max(y0, lb), ub);

        try
            [sol, info] = solve_trim_point_6dof_aero( ...
                Vt, alt_m, gamma_target_rad, db, y0, lb, ub, opts, reg);

            alpha_deg = sol(1);
            theta_rad = sol(2);
            PLA_deg   = sol(3);
            PFLAP_deg = sol(4);
            SPLRL_deg = sol(5);

            [thisFeasible, failCode, failReason] = classify_trim_result( ...
                info, PLA_deg, PFLAP_deg, SPLRL_deg, ...
                PLA_min, PLA_max, PFLAP_lim, SPLRL_lim);

            alpha_req(ia,iv) = alpha_deg;
            theta_req(ia,iv) = rad2deg(theta_rad);
            gamma_req(ia,iv) = info.gamma_deg;
            PLA_req(ia,iv)   = PLA_deg;
            PFLAP_req(ia,iv) = PFLAP_deg;
            SPLRL_req(ia,iv) = SPLRL_deg;
            feasible(ia,iv)  = thisFeasible;

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
                fprintf(['AERO-CLB | vy=%6.1f | alt=%6.0f | Vt=%6.1f | feas=%d | ' ...
                         'alpha=%7.3f | theta=%7.3f | gamma=%7.3f | ' ...
                         'PLA=%7.3f | PFLAP=%7.3f | SPLRL=%7.3f | res=%9.3e | exit=%d\n'], ...
                    vy, alt_m, Vt, thisFeasible, alpha_deg, rad2deg(theta_rad), ...
                    info.gamma_deg, PLA_deg, PFLAP_deg, SPLRL_deg, ...
                    info.resnorm_phys, info.exitflag);

                if ~thisFeasible
                    fprintf('           fail reason: %s\n', failReason);
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
                fprintf('AERO-CLB | vy=%6.1f | alt=%6.0f | Vt=%6.1f | feas=0 | ERROR: %s\n', ...
                    vy, alt_m, Vt, ME.message);
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
out.Vlist = Vlist;
out.alt_list_m = alt_list_m;
out.PLA_req = PLA_req;
out.PFLAP_req = PFLAP_req;
out.SPLRL_req = SPLRL_req;
out.alpha_req = alpha_req;
out.theta_req = theta_req;
out.gamma_req = gamma_req;
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
% CORE SOLVER
% =====================================================================

function [sol, info] = solve_trim_point_6dof_aero(Vt, alt_m, gamma_target_rad, db, y0, lb, ub, opts, reg)

obj = @(y) trim_residuals_6dof_safe_aero(y, Vt, alt_m, gamma_target_rad, db, reg);
[sol, resvec, ~, exitflag, output] = lsqnonlin(obj, y0, lb, ub, opts); %#ok<ASGLU>

[res_phys, x_trim, u_trim, extra] = trim_residuals_6dof_phys_aero(sol, Vt, alt_m, gamma_target_rad, db);
resnorm_phys = norm(res_phys);

du_tol = 5e-3;
dw_tol = 5e-3;
dq_tol = 5e-2;
dgamma_tol = deg2rad(0.3);

success = exitflag > 0 ...
    && all(isfinite(sol)) ...
    && all(isfinite(res_phys)) ...
    && abs(res_phys(1)) < du_tol ...
    && abs(res_phys(2)) < dw_tol ...
    && abs(res_phys(3)) < dq_tol ...
    && abs(res_phys(4)) < dgamma_tol;

info = struct();
info.exitflag = exitflag;
info.output = output;
info.residual = resvec;
info.residual_phys = res_phys;
info.resnorm = norm(resvec);
info.resnorm_phys = resnorm_phys;
info.success = success;
info.x_trim = x_trim;
info.u_trim = u_trim;
info.gamma_deg = rad2deg(extra.gamma_actual_rad);
info.gamma_err_deg = rad2deg(extra.gamma_actual_rad - gamma_target_rad);
info.dx = extra.dx;
end

function R = trim_residuals_6dof_safe_aero(y, Vt, alt_m, gamma_target_rad, db, reg)
try
    [Rphys, ~, ~, ~] = trim_residuals_6dof_phys_aero(y, Vt, alt_m, gamma_target_rad, db);
    if any(~isfinite(Rphys))
        R = 1e6 * ones(6,1);
        return
    end

    PFLAP_deg = y(4);
    SPLRL_deg = y(5);

    R = [ ...
        Rphys(:); ...
        sqrt(reg.w_pflap) * PFLAP_deg; ...
        sqrt(reg.w_splrl) * SPLRL_deg ...
    ];
catch
    R = 1e6 * ones(6,1);
end
end

function [R, x_trim, u_trim, extra] = trim_residuals_6dof_phys_aero(y, Vt, alt_m, gamma_target_rad, db)

alpha_deg = y(1);
theta_rad = y(2);
PLA_deg   = y(3);
PFLAP_deg = y(4);
SPLRL_deg = y(5);

if abs(alpha_deg) > 40 || ...
   abs(rad2deg(theta_rad)) > 50 || ...
   PLA_deg < -10 || PLA_deg > 150 || ...
   abs(PFLAP_deg) > 100 || ...
   abs(SPLRL_deg) > 100

    R = 1e6 * ones(4,1);
    x_trim = nan(12,1);
    u_trim = struct();
    extra = struct('gamma_actual_rad',NaN,'dx',nan(12,1));
    return
end

alpha_rad = deg2rad(alpha_deg);

u_b = Vt * cos(alpha_rad);
v_b = 0;
w_b = Vt * sin(alpha_rad);

x_trim = [u_b; v_b; w_b; 0; 0; 0; 0; theta_rad; 0; 0; 0; -alt_m];

u_trim = struct();
u_trim.PLA   = PLA_deg;
u_trim.PFLAP = PFLAP_deg;
u_trim.SPLRL = SPLRL_deg;

% Fixed-zero controls for this trim branch
u_trim.IBLFL = 0;
u_trim.OBLFL = 0;
u_trim.ELEVL = 0;
u_trim.CLAML = 0;
u_trim.AMTL  = 0;
u_trim.DEFL  = 0;

u_trim.IBLFR = 0;
u_trim.OBLFR = 0;
u_trim.ELEVR = 0;
u_trim.SPLRR = 0;
u_trim.CLAMR = 0;
u_trim.AMTR  = 0;
u_trim.DEFR  = 0;

u_trim.DELTAZ = 0;
u_trim.DELTAW = 0;

dx = ice_step_6dof_ctrl(0, x_trim, u_trim, db);
gamma_actual_rad = gamma_from_state12_local_aero(x_trim);

R = [dx(1); dx(3); dx(5); gamma_actual_rad - gamma_target_rad];

extra = struct();
extra.gamma_actual_rad = gamma_actual_rad;
extra.dx = dx;
end

% =====================================================================
% INITIAL GUESSES
% =====================================================================

function y0 = default_trim_guess_level_aero(Vt, alt_m, PLA_min, PLA_max)
alpha0_deg = max(1, min(10, 2 + 0.015*(200 - Vt) + 0.0001*alt_m));
theta0_rad = deg2rad(alpha0_deg);
PLA0_deg   = max(PLA_min + 5, min(PLA_max - 5, 50 + 0.004*alt_m + 0.04*(Vt - 100)));
PFLAP0_deg = 0;
SPLRL0_deg = 0;
y0 = [alpha0_deg; theta0_rad; PLA0_deg; PFLAP0_deg; SPLRL0_deg];
end

function y0 = default_trim_guess_climb_aero(Vt, alt_m, gamma_target_rad, PLA_min, PLA_max)
gamma_deg = rad2deg(gamma_target_rad);
alpha0_deg = max(1, min(12, 3 + 0.015*(200 - Vt) + 0.0001*alt_m));
theta0_rad = deg2rad(alpha0_deg + gamma_deg);
PLA0_deg   = max(PLA_min + 5, min(PLA_max - 5, 55 + 0.004*alt_m + 0.04*(Vt - 100) + 0.8*abs(gamma_deg)));
PFLAP0_deg = 0;
SPLRL0_deg = 0;
y0 = [alpha0_deg; theta0_rad; PLA0_deg; PFLAP0_deg; SPLRL0_deg];
end

% =====================================================================
% GAMMA
% =====================================================================

function gamma_rad = gamma_from_state12_local_aero(x12)
u = x12(1); v = x12(2); w = x12(3);
phi = x12(7); theta = x12(8); psi = x12(9);

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

function [isFeasible, code, reason] = classify_trim_result( ...
    info, PLA_deg, PFLAP_deg, SPLRL_deg, ...
    PLA_min, PLA_max, PFLAP_lim, SPLRL_lim)

% codes:
% 1 feasible
% 2 solver/lookup error
% 3 exitflag/non-convergence
% 4 PLA limit
% 5 PFLAP limit
% 6 SPLRL limit
% 7 du tolerance
% 8 dw tolerance
% 9 dq tolerance
% 10 gamma tolerance
% 11 other

du_tol = 5e-3;
dw_tol = 5e-3;
dq_tol = 5e-2;
dgamma_tol = deg2rad(0.3);

if info.success
    isFeasible = true;
    code = 1;
    reason = "Feasible";
    return
end

isFeasible = false;

if info.exitflag <= 0
    code = 3;
    reason = "Non-converged exitflag";
    return
end

if ~(PLA_deg >= PLA_min && PLA_deg <= PLA_max)
    code = 4;
    reason = "PLA limit";
    return
end

if abs(PFLAP_deg) > PFLAP_lim
    code = 5;
    reason = "PFLAP limit";
    return
end

if abs(SPLRL_deg) > SPLRL_lim
    code = 6;
    reason = "SPLRL limit";
    return
end

if abs(info.residual_phys(1)) >= du_tol
    code = 7;
    reason = "du tolerance";
    return
end

if abs(info.residual_phys(2)) >= dw_tol
    code = 8;
    reason = "dw tolerance";
    return
end

if abs(info.residual_phys(3)) >= dq_tol
    code = 9;
    reason = "dq tolerance";
    return
end

if abs(info.residual_phys(4)) >= dgamma_tol
    code = 10;
    reason = "gamma tolerance";
    return
end

code = 11;
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
    fprintf('\nNo border failures found for case "%s".\n', caseData.caseName);
    return
end

codes = caseData.failure_code(mask);
reasons = caseData.failure_reason(mask);

fprintf('\nBorder failure summary for case "%s":\n', caseData.caseName);

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
% CSV + PLOTTING
% =====================================================================

function save_case_csv_aero(caseData, stem)
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
    caseData.PFLAP_req(:), ...
    caseData.SPLRL_req(:), ...
    'VariableNames', { ...
    'Vt_mps','Altitude_m','Feasible','FailureCode','FailureReason', ...
    'Alpha_deg','Theta_deg','Gamma_deg', ...
    'PLA_deg','PFLAP_deg','SPLRL_deg'});

fname = matlab.lang.makeValidName(stem) + ".csv";
writetable(T, fname);
fprintf('Saved %s\n', fname);
end

function plot_case_aero(caseData, db, figTitle)
Vlist = caseData.Vlist;
alt_list_m = caseData.alt_list_m;
PLA_req   = caseData.PLA_req;
PFLAP_req = caseData.PFLAP_req;
SPLRL_req = caseData.SPLRL_req;
feasible  = caseData.feasible;
failure_code = caseData.failure_code;
border_fail_mask = caseData.border_fail_mask;

[Vgrid, AltGrid] = meshgrid(Vlist, alt_list_m);

PLA_min = min(db.eng.PLA);
PLA_max = max(db.eng.PLA);
PFLAP_lim = max(abs(db.aero.ctrl.PFLAP));
SPLRL_lim = max(abs(db.aero.ctrl.SPLRL2));

canSurf = (size(Vgrid,1) >= 2) && (size(Vgrid,2) >= 2);

if canSurf
    plot_surface_case(Vgrid, AltGrid, PLA_req, feasible, PLA_min, PLA_max, ...
        'PLA (deg)', sprintf('%s: Required PLA', figTitle), 'PLA');
    plot_surface_case(Vgrid, AltGrid, PFLAP_req, feasible, -PFLAP_lim, PFLAP_lim, ...
        'PFLAP (deg)', sprintf('%s: Required PFLAP', figTitle), 'PFLAP');
    plot_surface_case(Vgrid, AltGrid, SPLRL_req, feasible, -SPLRL_lim, SPLRL_lim, ...
        'SPLRL (deg)', sprintf('%s: Required SPLRL', figTitle), 'SPLRL');
else
    plot_line_case(Vlist, alt_list_m, PLA_req, feasible, 'PLA (deg)', sprintf('%s: Required PLA', figTitle));
    plot_line_case(Vlist, alt_list_m, PFLAP_req, feasible, 'PFLAP (deg)', sprintf('%s: Required PFLAP', figTitle));
    plot_line_case(Vlist, alt_list_m, SPLRL_req, feasible, 'SPLRL (deg)', sprintf('%s: Required SPLRL', figTitle));
end

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

display_code = zeros(size(failure_code));
display_code(feasible) = 1;
display_code(border_fail_mask) = failure_code(border_fail_mask);

figure('Color','w');
imagesc(Vlist, alt_list_m, display_code);
set(gca,'YDir','normal')

cmap = [ ...
    1.00 1.00 1.00;   % 0 not shown / non-border infeasible
    0.00 0.60 0.00;   % 1 feasible
    0.00 0.00 0.00;   % 2 solver
    0.45 0.45 0.45;   % 3 exitflag
    0.85 0.10 0.10;   % 4 PLA
    1.00 0.50 0.00;   % 5 PFLAP
    1.00 0.85 0.10;   % 6 SPLRL
    0.10 0.30 0.90;   % 7 du
    0.00 0.75 0.90;   % 8 dw
    0.85 0.00 0.85;   % 9 dq
    0.55 0.20 0.75;   % 10 gamma
    0.55 0.35 0.15];  % 11 other

colormap(gca, cmap)
caxis([0 11])
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: Border Failure Reasons', figTitle))

cb = colorbar;
cb.Ticks = 0:11;
cb.TickLabels = { ...
    'Not shown', ...
    'Feasible', ...
    'Solver error', ...
    'Exitflag', ...
    'PLA limit', ...
    'PFLAP limit', ...
    'SPLRL limit', ...
    'du tol', ...
    'dw tol', ...
    'dq tol', ...
    'gamma tol', ...
    'Other'};
end

function plot_surface_case(Vgrid, AltGrid, Zreq, feasible, zlo, zhi, zlabeltxt, titletxt, legbase)
figure('Color','w');
hold on

Z_feas = Zreq;
Z_infeas = Zreq;
Z_feas(~feasible) = NaN;
Z_infeas(feasible) = NaN;

s1 = surf(Vgrid, AltGrid, Z_feas, Z_feas, 'EdgeColor','none', 'FaceColor','interp');

colormap(gca, turbo)
cb = colorbar;
cb.Label.String = zlabeltxt;

viol = zeros(size(Zreq));
viol(Zreq > zhi) = Zreq(Zreq > zhi) - zhi;
viol(Zreq < zlo) = zlo - Zreq(Zreq < zlo);
viol(feasible) = NaN;

redmap = zeros([size(Zreq), 3]);
if any(~isnan(viol(:)))
    vmax = max(viol(:), [], 'omitnan');
    if vmax == 0, vmax = 1; end
    vnorm = viol / vmax;
    redmap(:,:,1) = 1.0;
    redmap(:,:,2) = 0.9 - 0.7*vnorm;
    redmap(:,:,3) = 0.9 - 0.7*vnorm;
else
    redmap(:,:,1) = 1.0;
    redmap(:,:,2) = 0.9;
    redmap(:,:,3) = 0.9;
end

s2 = surf(Vgrid, AltGrid, Z_infeas, redmap, 'EdgeColor','none', 'FaceColor','interp');
s3 = surf(Vgrid, AltGrid, zhi*ones(size(Zreq)), 'FaceAlpha',0.18, 'EdgeColor','none', 'FaceColor',[0 0 0]);
s4 = surf(Vgrid, AltGrid, zlo*ones(size(Zreq)), 'FaceAlpha',0.18, 'EdgeColor','none', 'FaceColor',[0 0 0]);

grid on
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
zlabel(zlabeltxt)
title(titletxt)
view(45,30)

legend([s1 s2 s3 s4], ...
    {['Feasible ' legbase], ['Out-of-bounds ' legbase], [legbase ' upper limit'], [legbase ' lower limit']}, ...
    'Location','best');
end

function plot_line_case(Vlist, alt_list_m, Zreq, feasible, ylab, titletxt)
figure('Color','w');
hold on

z = Zreq(:);
f = feasible(:);

if numel(Vlist) == numel(z)
    x = Vlist(:);
    xlabeltxt = 'V_t (m/s)';
elseif numel(alt_list_m) == numel(z)
    x = alt_list_m(:);
    xlabeltxt = 'Altitude (m)';
else
    x = (1:numel(z)).';
    xlabeltxt = 'Index';
end

plot(x(f), z(f), 'o-', 'LineWidth', 1.4);
plot(x(~f), z(~f), 'x', 'LineWidth', 1.2);
grid on
xlabel(xlabeltxt)
ylabel(ylab)
title(titletxt)
end