function results = sweep_trim_vt_6dof_aero_equilibrium_map(mode)
% SWEEP_TRIM_VT_6DOF_AERO_EQUILIBRIUM_MAP
%
% Aero-control 6DOF equilibrium sweep written in the same overall style as
% sweep_trim_vt_6dof_equilibrium_map.m, but using the ORIGINAL aero trim
% acceptance logic.
%
% Trim variables:
%   y = [alpha_deg; theta_rad; PLA_deg; PFLAP_deg; SPLRL_deg]
%
% Important:
% - Equilibrium acceptance uses the original aero component tolerances:
%       |du|     < 5e-3
%       |dw|     < 5e-3
%       |dq|     < 5e-2
%       |dgamma| < 0.3 deg
% - No extra actuator-limit feasibility logic is applied in the output.
% - Numerical bounds are still used for solver robustness.
%
% Output CSV is intended to support the same PLA contour workflow as the
% TVC equilibrium-map code.
%
% Modes:
%   'level' only

if nargin < 1
    mode = 'level';
end

mode = lower(string(mode));
db = ice_init_aeroctrl();

% ---------------- USER OPTIONS ----------------
Vlist = 480:4:560;               % m/s
alt_list_m = 0:200:20000;     % m

showConsoleProgress = true;
saveCSV = true;
makePlots = true;

nConsecutiveNoEqStop = 1;
pointTimeout_s = 1.0;

% Regularisation copied from original aero code
reg = struct();
reg.w_pflap = 1e-3;
reg.w_splrl = 1e-3;

% Numerical bounds only (not reported as feasibility limits)
alpha_min_deg = -40;
alpha_max_deg =  40;

theta_min_deg = -50;
theta_max_deg =  50;

PLA_min_num = min(db.eng.PLA);
PLA_max_num = max(db.eng.PLA) + 1000;

PFLAP_lim_num = max(abs(db.aero.ctrl.PFLAP));
SPLRL_lim_num = max(abs(db.aero.ctrl.SPLRL2));
% ---------------------------------------------

results = struct();
results.settings = struct( ...
    'Vlist', Vlist, ...
    'alt_list_m', alt_list_m, ...
    'nConsecutiveNoEqStop', nConsecutiveNoEqStop, ...
    'pointTimeout_s', pointTimeout_s, ...
    'alpha_min_deg', alpha_min_deg, ...
    'alpha_max_deg', alpha_max_deg, ...
    'theta_min_deg', theta_min_deg, ...
    'theta_max_deg', theta_max_deg, ...
    'PLA_min_num', PLA_min_num, ...
    'PLA_max_num', PLA_max_num, ...
    'PFLAP_lim_num', PFLAP_lim_num, ...
    'SPLRL_lim_num', SPLRL_lim_num);

if mode ~= "level"
    error('This version currently supports only mode = "level".');
end

results.level = run_level_case_6dof_aero_equilibrium_map( ...
    db, Vlist, alt_list_m, ...
    alpha_min_deg, alpha_max_deg, ...
    theta_min_deg, theta_max_deg, ...
    PLA_min_num, PLA_max_num, ...
    PFLAP_lim_num, SPLRL_lim_num, ...
    showConsoleProgress, nConsecutiveNoEqStop, pointTimeout_s, reg);

if makePlots
    plot_case_aero_equilibrium_map(results.level, ...
        '6-DOF Aero-Control Equilibrium Map');
end

if saveCSV
    save_case_csv_aero_equilibrium_map(results.level, 'level_6dof_aero_equilibrium_map_extra');
end

end

% =====================================================================
% LEVEL CASE
% =====================================================================

function out = run_level_case_6dof_aero_equilibrium_map( ...
    db, Vlist, alt_list_m, ...
    alpha_min_deg, alpha_max_deg, ...
    theta_min_deg, theta_max_deg, ...
    PLA_min_num, PLA_max_num, ...
    PFLAP_lim_num, SPLRL_lim_num, ...
    showConsoleProgress, nConsecutiveNoEqStop, pointTimeout_s, reg)

nV = numel(Vlist);
nA = numel(alt_list_m);

PLA_req   = NaN(nA,nV);
PFLAP_req = NaN(nA,nV);
SPLRL_req = NaN(nA,nV);
alpha_req = NaN(nA,nV);
theta_req = NaN(nA,nV);
gamma_req = NaN(nA,nV);

equilibrium_ok = false(nA,nV);

failure_code   = zeros(nA,nV);
failure_reason = strings(nA,nV);
failure_reason(:) = "Not evaluated";
evaluated      = false(nA,nV);

exitflag_map = NaN(nA,nV);
resnorm_map  = NaN(nA,nV);

x_trim_cell = cell(nA,nV);
u_trim_cell = cell(nA,nV);
trim_info   = cell(nA,nV);

opts_base = optimoptions('lsqnonlin', ...
    'Display','off', ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10, ...
    'MaxFunctionEvaluations',2500, ...
    'MaxIterations',300);

lb = [ ...
    alpha_min_deg; ...
    deg2rad(theta_min_deg); ...
    PLA_min_num; ...
   -PFLAP_lim_num; ...
   -SPLRL_lim_num];

ub = [ ...
    alpha_max_deg; ...
    deg2rad(theta_max_deg); ...
    PLA_max_num; ...
    PFLAP_lim_num; ...
    SPLRL_lim_num];

for iv = 1:nV
    Vt = Vlist(iv);
    prev_sol = [];

    sawEquilibrium = false;
    consecutiveNoEq = 0;

    for ia = 1:nA
        alt_m = alt_list_m(ia);

        if isempty(prev_sol)
            y0 = default_trim_guess_level_aero_equilibrium(Vt, alt_m, PLA_min_num, PLA_max_num);
        else
            y0 = prev_sol;
        end
        y0 = min(max(y0, lb), ub);

        solveStart = tic;
        opts = optimoptions(opts_base, ...
            'OutputFcn', @(x,optimValues,state) local_timeout_outputfcn(x,optimValues,state,solveStart,pointTimeout_s));

        try
            [sol, info] = solve_trim_point_6dof_aero_equilibrium( ...
                Vt, alt_m, 0.0, db, y0, lb, ub, opts, reg, pointTimeout_s);

            alpha_deg = sol(1);
            theta_rad = sol(2);
            PLA_deg   = sol(3);
            PFLAP_deg = sol(4);
            SPLRL_deg = sol(5);

            [thisEq, failCode, failReason] = classify_trim_result_equilibrium_only_aero(info, pointTimeout_s);

            alpha_req(ia,iv) = alpha_deg;
            theta_req(ia,iv) = rad2deg(theta_rad);
            gamma_req(ia,iv) = info.gamma_deg;
            PLA_req(ia,iv)   = PLA_deg;
            PFLAP_req(ia,iv) = PFLAP_deg;
            SPLRL_req(ia,iv) = SPLRL_deg;

            equilibrium_ok(ia,iv) = thisEq;
            evaluated(ia,iv)      = true;
            failure_code(ia,iv)   = failCode;
            failure_reason(ia,iv) = failReason;
            exitflag_map(ia,iv)   = info.exitflag;
            resnorm_map(ia,iv)    = info.resnorm_phys;

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
                fprintf(['AERO-EQ | alt=%6.0f | Vt=%6.1f | eq=%d | ' ...
                         'alpha=%7.3f | theta=%7.3f | gamma=%7.3f | ' ...
                         'PLA=%9.3f | PFLAP=%8.3f | SPLRL=%8.3f | ' ...
                         'res=%9.3e | exit=%d | t=%6.2f s | %s\n'], ...
                    alt_m, Vt, thisEq, ...
                    alpha_deg, rad2deg(theta_rad), info.gamma_deg, ...
                    PLA_deg, PFLAP_deg, SPLRL_deg, ...
                    info.resnorm_phys, info.exitflag, info.solveTime_s, failReason);
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
                fprintf('AERO-EQ | alt=%6.0f | Vt=%6.1f | eq=0 | ERROR: %s\n', ...
                    alt_m, Vt, ME.message);
            end
        end

        if sawEquilibrium && consecutiveNoEq >= nConsecutiveNoEqStop
            if showConsoleProgress
                fprintf(['AERO-EQ | alt=%6.0f | Vt=%6.1f | STOP after %d consecutive no-equilibrium points\n'], ...
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

out.PLA_req   = PLA_req;
out.PFLAP_req = PFLAP_req;
out.SPLRL_req = SPLRL_req;
out.alpha_req = alpha_req;
out.theta_req = theta_req;
out.gamma_req = gamma_req;

out.equilibrium_ok = equilibrium_ok;
out.failure_code = failure_code;
out.failure_reason = failure_reason;
out.evaluated = evaluated;
out.exitflag_map = exitflag_map;
out.resnorm_map = resnorm_map;

out.x_trim_cell = x_trim_cell;
out.u_trim_cell = u_trim_cell;
out.trim_info = trim_info;

end

% =====================================================================
% CORE SOLVER
% =====================================================================

function [sol, info] = solve_trim_point_6dof_aero_equilibrium( ...
    Vt, alt_m, gamma_target_rad, db, y0, lb, ub, opts, reg, pointTimeout_s)

obj = @(y) trim_residuals_6dof_safe_aero_equilibrium(y, Vt, alt_m, gamma_target_rad, db, reg);

tSolve = tic;
[sol, resvec, ~, exitflag, output] = lsqnonlin(obj, y0, lb, ub, opts); %#ok<ASGLU>
solveTime_s = toc(tSolve);

[res_phys, x_trim, u_trim, extra] = trim_residuals_6dof_phys_aero_equilibrium(sol, Vt, alt_m, gamma_target_rad, db);
resnorm_phys = norm(res_phys);

du_tol = 5e-3;
dw_tol = 5e-3;
dq_tol = 5e-2;
dgamma_tol = deg2rad(0.3);

timedOut = (solveTime_s >= pointTimeout_s - 1e-3) && (exitflag <= 0);

success = ~timedOut ...
    && exitflag > 0 ...
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
info.timedOut = timedOut;
info.solveTime_s = solveTime_s;
info.x_trim = x_trim;
info.u_trim = u_trim;
info.gamma_deg = rad2deg(extra.gamma_actual_rad);
info.gamma_err_deg = rad2deg(extra.gamma_actual_rad - gamma_target_rad);
info.dx = extra.dx;

end

function R = trim_residuals_6dof_safe_aero_equilibrium(y, Vt, alt_m, gamma_target_rad, db, reg)
try
    [Rphys, ~, ~, ~] = trim_residuals_6dof_phys_aero_equilibrium(y, Vt, alt_m, gamma_target_rad, db);

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

function [R, x_trim, u_trim, extra] = trim_residuals_6dof_phys_aero_equilibrium(y, Vt, alt_m, gamma_target_rad, db)

% y = [alpha_deg; theta_rad; PLA_deg; PFLAP_deg; SPLRL_deg]

alpha_deg = y(1);
theta_rad = y(2);
PLA_deg   = y(3);
PFLAP_deg = y(4);
SPLRL_deg = y(5);

if abs(alpha_deg) > 40 || ...
   abs(rad2deg(theta_rad)) > 50 || ...
   PLA_deg < -20 || PLA_deg > 1000 || ...
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
gamma_actual_rad = gamma_from_state12_local_aero_equilibrium(x_trim);

R = [dx(1); dx(3); dx(5); gamma_actual_rad - gamma_target_rad];

extra = struct();
extra.gamma_actual_rad = gamma_actual_rad;
extra.dx = dx;

end

% =====================================================================
% INITIAL GUESS
% =====================================================================

function y0 = default_trim_guess_level_aero_equilibrium(Vt, alt_m, PLA_min_num, PLA_max_num)

alpha0_deg = max(1, min(10, 2 + 0.015*(200 - Vt) + 0.0001*alt_m));
theta0_rad = deg2rad(alpha0_deg);
PLA0_deg   = min(max(50 + 0.004*alt_m + 0.04*(Vt - 100), PLA_min_num + 5), PLA_max_num - 5);
PFLAP0_deg = 0;
SPLRL0_deg = 0;

y0 = [alpha0_deg; theta0_rad; PLA0_deg; PFLAP0_deg; SPLRL0_deg];

end

% =====================================================================
% GAMMA
% =====================================================================

function gamma_rad = gamma_from_state12_local_aero_equilibrium(x12)

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
% EQUILIBRIUM CLASSIFICATION
% =====================================================================

function [isEquilibrium, code, reason] = classify_trim_result_equilibrium_only_aero(info, pointTimeout_s)

du_tol = 5e-3;
dw_tol = 5e-3;
dq_tol = 5e-2;
dgamma_tol = deg2rad(0.3);

isEquilibrium = false;

if isfield(info,'timedOut') && info.timedOut
    code = 9;
    reason = sprintf('Timeout (>%g s)', pointTimeout_s);
    return
end

if info.exitflag <= 0
    code = 3;
    reason = "Non-converged exitflag";
    return
end

if ~all(isfinite(info.residual_phys))
    code = 4;
    reason = "Non-finite residual";
    return
end

if abs(info.residual_phys(1)) >= du_tol
    code = 5;
    reason = "du tolerance";
    return
end

if abs(info.residual_phys(2)) >= dw_tol
    code = 6;
    reason = "dw tolerance";
    return
end

if abs(info.residual_phys(3)) >= dq_tol
    code = 7;
    reason = "dq tolerance";
    return
end

if abs(info.residual_phys(4)) >= dgamma_tol
    code = 8;
    reason = "gamma tolerance";
    return
end

isEquilibrium = true;
code = 1;
reason = "Equilibrium";

end

% =====================================================================
% CSV EXPORT
% =====================================================================

function save_case_csv_aero_equilibrium_map(caseData, stem)

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
    caseData.PFLAP_req(:), ...
    caseData.SPLRL_req(:), ...
    'VariableNames', { ...
    'Vt_mps','Altitude_m','EquilibriumOK','FailureCode','FailureReason', ...
    'Exitflag','Resnorm', ...
    'Alpha_deg','Theta_deg','Gamma_deg', ...
    'PLA_deg','PFLAP_deg','SPLRL_deg'});

fname = matlab.lang.makeValidName(stem) + ".csv";
writetable(T, fname);
fprintf('Saved %s\n', fname);

end

% =====================================================================
% PLOTTING
% =====================================================================

function plot_case_aero_equilibrium_map(caseData, figTitle)

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
title(sprintf('%s: PLA', figTitle))
colorbar
colormap(gca, turbo)

figure('Color','w');
imagesc(Vlist, alt_list_m, caseData.PFLAP_req);
set(gca,'YDir','normal')
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: PFLAP', figTitle))
colorbar
colormap(gca, turbo)

figure('Color','w');
imagesc(Vlist, alt_list_m, caseData.SPLRL_req);
set(gca,'YDir','normal')
xlabel('V_t (m/s)')
ylabel('Altitude (m)')
title(sprintf('%s: SPLRL', figTitle))
colorbar
colormap(gca, turbo)

end

% =====================================================================
% TIMEOUT OUTPUT FUNCTION
% =====================================================================

function stop = local_timeout_outputfcn(~, ~, state, solveStart, pointTimeout_s)
stop = false;

switch char(state)
    case {'init','iter','interrupt'}
        if toc(solveStart) > pointTimeout_s
            stop = true;
        end
end

end