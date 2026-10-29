function results = tvc_turn_envelope_6dof_tas()
% TVC_TURN_ENVELOPE_6DOF_PSIDOT_TAS
%
% Representative 6DOF coordinated-turn trim envelope for TVC-only control,
% built by searching maximum steady heading rate at each (Vt, altitude).
%
% At each (Vt, h), the code:
%   1) searches over candidate steady heading rate psi_dot
%   2) solves a turning trim consistent with that psi_dot
%   3) records the largest feasible trimmed turn
%   4) converts that to n, phi, turn radius
%
% Unknowns at each fixed (Vt,h,psi_dot):
%   y = [alpha_deg, beta_deg, theta_rad, phi_rad, p_rad_s, q_rad_s, r_rad_s, ...
%        PLA_deg, DELTAZ_deg, DELTAW_deg]
%
% Residuals:
%   du = 0
%   dv = 0
%   dw = 0
%   dp = 0
%   dq = 0
%   dr = 0
%   gamma = 0
%   phi_dot = 0
%   theta_dot = 0
%   psi_dot - psi_dot_target = 0
%
% Notes:
% - This is a practical steady-turn trim approximation using the full 6DOF model.
% - It is much more meaningful than prescribing bank angle alone.
% - The resulting V-n graph is representative rather than a formal EM diagram.

db = ice_init();

% ============================================================
% USER SETTINGS
% ============================================================
alt_list_m = [0, 10000, 20000];

vt_list_by_alt = {
    76:2:330   % for 0 m
    210:2:460   % for 10000 m
    460:2:550   % for 20000 m
};

PLA_min = min(db.eng.PLA);
PLA_max = max(db.eng.PLA);

DELTAZ_lim = 25;      % deg
DELTAW_lim = 25;      % deg
alpha_lim  = 20;      % deg
beta_lim   = 8;       % deg
theta_lim_deg = 35;   % deg
phi_lim_deg   = 85;   % deg
p_lim_deg_s   = 60;   % deg/s
q_lim_deg_s   = 60;   % deg/s
r_lim_deg_s   = 120;  % deg/s

psidot_min_deg_s = 0;
psidot_max_deg_s = 25;
psidot_step_deg_s = 1.0;
psidot_tol_deg_s  = 0.1;

showConsoleProgress = true;
saveCSV = true;
makePlots = true;

trim_tol = 1e-5;
vt_point_time_limit_s = 5.0;   % max total time spent on one (alt,Vt) point
% ============================================================

nAlt = numel(alt_list_m);

if numel(vt_list_by_alt) ~= nAlt
    error('vt_list_by_alt must have one entry per altitude in alt_list_m.');
end

nVt_max = max(cellfun(@numel, vt_list_by_alt));

Vt_grid_mps       = NaN(nAlt, nVt_max);
psidot_max_map    = NaN(nAlt, nVt_max);
phi_deg_map       = NaN(nAlt, nVt_max);
n_max_map         = NaN(nAlt, nVt_max);
radius_m_map      = NaN(nAlt, nVt_max);

alpha_deg_map     = NaN(nAlt, nVt_max);
beta_deg_map      = NaN(nAlt, nVt_max);
theta_deg_map     = NaN(nAlt, nVt_max);
p_deg_s_map       = NaN(nAlt, nVt_max);
q_deg_s_map       = NaN(nAlt, nVt_max);
r_deg_s_map       = NaN(nAlt, nVt_max);
PLA_deg_map       = NaN(nAlt, nVt_max);
DELTAZ_deg_map    = NaN(nAlt, nVt_max);
DELTAW_deg_map    = NaN(nAlt, nVt_max);
TAS_map           = NaN(nAlt, nVt_max);
resnorm_map       = NaN(nAlt, nVt_max);

feasible_map      = false(nAlt, nVt_max);

failure_code_map   = zeros(nAlt, nVt_max);
failure_reason_map = strings(nAlt, nVt_max);
failure_reason_map(:) = "Not evaluated";

psidot_fail_map    = NaN(nAlt, nVt_max);
alpha_fail_map     = NaN(nAlt, nVt_max);
beta_fail_map      = NaN(nAlt, nVt_max);
theta_fail_map     = NaN(nAlt, nVt_max);
phi_fail_map       = NaN(nAlt, nVt_max);
p_fail_map         = NaN(nAlt, nVt_max);
q_fail_map         = NaN(nAlt, nVt_max);
r_fail_map         = NaN(nAlt, nVt_max);
PLA_fail_map       = NaN(nAlt, nVt_max);
DELTAZ_fail_map    = NaN(nAlt, nVt_max);
DELTAW_fail_map    = NaN(nAlt, nVt_max);
resnorm_fail_map   = NaN(nAlt, nVt_max);

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-9, ...
    'StepTolerance',1e-9, ...
    'MaxFunctionEvaluations',3000, ...
    'MaxIterations',300);

for ia = 1:nAlt
    alt_m = alt_list_m(ia);
    vt_list = vt_list_by_alt{ia};
    nVt = numel(vt_list);

    yguess_speed = [];

    for iv = 1:nVt
        Vt = vt_list(iv);
        Vt_grid_mps(ia,iv) = Vt;
        TAS_map(ia,iv) = Vt;

        if isempty(yguess_speed)
            yguess_base = default_turn_guess(Vt, alt_m, PLA_min, PLA_max);
        else
            yguess_base = yguess_speed;
        end

        [psidot_best_deg_s, y_best, info_best, ok, failInfo] = find_max_psidot_at_point( ...
            Vt, alt_m, db, opts, trim_tol, ...
            yguess_base, ...
            PLA_min, PLA_max, DELTAZ_lim, DELTAW_lim, ...
            alpha_lim, beta_lim, theta_lim_deg, phi_lim_deg, ...
            p_lim_deg_s, q_lim_deg_s, r_lim_deg_s, ...
            psidot_min_deg_s, psidot_max_deg_s, psidot_step_deg_s, psidot_tol_deg_s, ...
            vt_point_time_limit_s);

        failure_code_map(ia,iv)   = failInfo.code;
        failure_reason_map(ia,iv) = failInfo.reason;

        psidot_fail_map(ia,iv)    = failInfo.psidot_fail_deg_s;
        alpha_fail_map(ia,iv)     = failInfo.y_fail(1);
        beta_fail_map(ia,iv)      = failInfo.y_fail(2);
        theta_fail_map(ia,iv)     = rad2deg_safe(failInfo.y_fail(3));
        phi_fail_map(ia,iv)       = rad2deg_safe(failInfo.y_fail(4));
        p_fail_map(ia,iv)         = rad2deg_safe(failInfo.y_fail(5));
        q_fail_map(ia,iv)         = rad2deg_safe(failInfo.y_fail(6));
        r_fail_map(ia,iv)         = rad2deg_safe(failInfo.y_fail(7));
        PLA_fail_map(ia,iv)       = failInfo.y_fail(8);
        DELTAZ_fail_map(ia,iv)    = failInfo.y_fail(9);
        DELTAW_fail_map(ia,iv)    = failInfo.y_fail(10);
        resnorm_fail_map(ia,iv)   = failInfo.resnorm_fail;

        if ok
            phi_rad = y_best(4);
            psidot_rad_s = deg2rad(psidot_best_deg_s);

            n_val = 1 / max(cos(phi_rad), 1e-9);

            if abs(psidot_rad_s) > 1e-9
                radius_m = Vt / abs(psidot_rad_s);
            else
                radius_m = Inf;
            end

            psidot_max_map(ia,iv)  = psidot_best_deg_s;
            phi_deg_map(ia,iv)     = rad2deg(phi_rad);
            n_max_map(ia,iv)       = n_val;
            radius_m_map(ia,iv)    = radius_m;

            alpha_deg_map(ia,iv)   = y_best(1);
            beta_deg_map(ia,iv)    = y_best(2);
            theta_deg_map(ia,iv)   = rad2deg(y_best(3));
            p_deg_s_map(ia,iv)     = rad2deg(y_best(5));
            q_deg_s_map(ia,iv)     = rad2deg(y_best(6));
            r_deg_s_map(ia,iv)     = rad2deg(y_best(7));
            PLA_deg_map(ia,iv)     = y_best(8);
            DELTAZ_deg_map(ia,iv)  = y_best(9);
            DELTAW_deg_map(ia,iv)  = y_best(10);
            resnorm_map(ia,iv)     = info_best.resnorm;

            feasible_map(ia,iv) = true;
            yguess_speed = y_best;

            if showConsoleProgress
                fprintf(['alt=%6.0f | TAS=%7.2f | feas=1 | ' ...
                         'psidot=%6.2f deg/s | n=%6.3f | phi=%6.2f | R=%9.1f | ' ...
                         'alpha=%6.2f | beta=%6.2f | theta=%6.2f | ' ...
                         'p=%6.2f | q=%6.2f | r=%6.2f | ' ...
                         'PLA=%7.2f | DZ=%6.2f | DW=%6.2f | res=%9.2e | fail=%s\n'], ...
                    alt_m, Vt, psidot_best_deg_s, n_val, rad2deg(phi_rad), radius_m, ...
                    y_best(1), y_best(2), rad2deg(y_best(3)), ...
                    rad2deg(y_best(5)), rad2deg(y_best(6)), rad2deg(y_best(7)), ...
                    y_best(8), y_best(9), y_best(10), info_best.resnorm, failInfo.reason);
            end
        else
            if showConsoleProgress
                fprintf('alt=%6.0f | TAS=%7.2f | feas=0 | fail=%s\n', ...
                    alt_m, Vt, failInfo.reason);
            end
        end
    end
end

results = struct();
results.altitude_m            = alt_list_m(:);
results.vt_list_by_alt        = vt_list_by_alt;
results.Vt_grid_mps           = Vt_grid_mps;
results.TAS_map               = TAS_map;
results.psidot_max_deg_s_map  = psidot_max_map;
results.phi_deg_map           = phi_deg_map;
results.n_max_map             = n_max_map;
results.radius_m_map          = radius_m_map;

results.alpha_deg_map    = alpha_deg_map;
results.beta_deg_map     = beta_deg_map;
results.theta_deg_map    = theta_deg_map;
results.p_deg_s_map      = p_deg_s_map;
results.q_deg_s_map      = q_deg_s_map;
results.r_deg_s_map      = r_deg_s_map;
results.PLA_deg_map      = PLA_deg_map;
results.DELTAZ_deg_map   = DELTAZ_deg_map;
results.DELTAW_deg_map   = DELTAW_deg_map;
results.resnorm_map      = resnorm_map;
results.feasible_map     = feasible_map;

results.failure_code_map   = failure_code_map;
results.failure_reason_map = failure_reason_map;
results.psidot_fail_map    = psidot_fail_map;
results.alpha_fail_map     = alpha_fail_map;
results.beta_fail_map      = beta_fail_map;
results.theta_fail_map     = theta_fail_map;
results.phi_fail_map       = phi_fail_map;
results.p_fail_map         = p_fail_map;
results.q_fail_map         = q_fail_map;
results.r_fail_map         = r_fail_map;
results.PLA_fail_map       = PLA_fail_map;
results.DELTAZ_fail_map    = DELTAZ_fail_map;
results.DELTAW_fail_map    = DELTAW_fail_map;
results.resnorm_fail_map   = resnorm_fail_map;

% ============================================================
% PLOTS
% ============================================================
if makePlots
    figure('Color','w'); hold on
    for ia = 1:nAlt
        y = n_max_map(ia,:);
        y(~feasible_map(ia,:)) = NaN;
        vt_row = Vt_grid_mps(ia,:);
        mask = isfinite(vt_row) & isfinite(y);
        plot(vt_row(mask), y(mask), 'LineWidth', 2, ...
            'DisplayName', sprintf('%g m', alt_list_m(ia)));
    end
    grid on
    xlabel('True airspeed V_t (m/s)')
    ylabel('Maximum trimmed load factor n')
    title('TVC-only representative V-n graph')
    legend('Location','best')

    figure('Color','w'); hold on
    for ia = 1:nAlt
        y = psidot_max_map(ia,:);
        y(~feasible_map(ia,:)) = NaN;
        vt_row = Vt_grid_mps(ia,:);
        mask = isfinite(vt_row) & isfinite(y);
        plot(vt_row(mask), y(mask), 'LineWidth', 2, ...
            'DisplayName', sprintf('%g m', alt_list_m(ia)));
    end
    grid on
    xlabel('True airspeed V_t (m/s)')
    ylabel('Maximum trimmed heading rate \psi̇ (deg/s)')
    title('TVC-only representative turn-rate graph')
    legend('Location','best')

    figure('Color','w'); hold on
    for ia = 1:nAlt
        y = phi_deg_map(ia,:);
        y(~feasible_map(ia,:)) = NaN;
        vt_row = Vt_grid_mps(ia,:);
        mask = isfinite(vt_row) & isfinite(y);
        plot(vt_row(mask), y(mask), 'LineWidth', 2, ...
            'DisplayName', sprintf('%g m', alt_list_m(ia)));
    end
    grid on
    xlabel('True airspeed V_t (m/s)')
    ylabel('Trim bank angle \phi (deg)')
    title('TVC-only representative bank-angle graph')
    legend('Location','best')

    figure('Color','w');
    tiledlayout(3,1)

    nexttile; hold on
    for ia = 1:nAlt
        y = PLA_deg_map(ia,:);
        y(~feasible_map(ia,:)) = NaN;
        vt_row = Vt_grid_mps(ia,:);
        mask = isfinite(vt_row) & isfinite(y);
        plot(vt_row(mask), y(mask), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('%g m', alt_list_m(ia)));
    end
    yline(PLA_min, '--k', 'PLA min');
    yline(PLA_max, '--k', 'PLA max');
    grid on
    xlabel('True airspeed V_t (m/s)')
    ylabel('PLA (deg)')
    title('Witness PLA')

    nexttile; hold on
    for ia = 1:nAlt
        y = DELTAZ_deg_map(ia,:);
        y(~feasible_map(ia,:)) = NaN;
        vt_row = Vt_grid_mps(ia,:);
        mask = isfinite(vt_row) & isfinite(y);
        plot(vt_row(mask), y(mask), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('%g m', alt_list_m(ia)));
    end
    yline(-DELTAZ_lim, '--k');
    yline( DELTAZ_lim, '--k');
    grid on
    xlabel('True airspeed V_t (m/s)')
    ylabel('\Delta Z (deg)')
    title('Witness \Delta Z')

    nexttile; hold on
    for ia = 1:nAlt
        y = DELTAW_deg_map(ia,:);
        y(~feasible_map(ia,:)) = NaN;
        vt_row = Vt_grid_mps(ia,:);
        mask = isfinite(vt_row) & isfinite(y);
        plot(vt_row(mask), y(mask), 'LineWidth', 1.6, ...
            'DisplayName', sprintf('%g m', alt_list_m(ia)));
    end
    yline(-DELTAW_lim, '--k');
    yline( DELTAW_lim, '--k');
    grid on
    xlabel('True airspeed V_t (m/s)')
    ylabel('\Delta W (deg)')
    title('Witness \Delta W')
    legend('Location','best')
end

% ============================================================
% CSV
% ============================================================
if saveCSV
    ALT_grid = repmat(alt_list_m(:), 1, nVt_max);
    VT_grid = Vt_grid_mps;

    T = table( ...
        ALT_grid(:), VT_grid(:), TAS_map(:), feasible_map(:), ...
        psidot_max_map(:), phi_deg_map(:), n_max_map(:), radius_m_map(:), ...
        alpha_deg_map(:), beta_deg_map(:), theta_deg_map(:), ...
        p_deg_s_map(:), q_deg_s_map(:), r_deg_s_map(:), ...
        PLA_deg_map(:), DELTAZ_deg_map(:), DELTAW_deg_map(:), ...
        resnorm_map(:), ...
        failure_code_map(:), cellstr(failure_reason_map(:)), ...
        psidot_fail_map(:), alpha_fail_map(:), beta_fail_map(:), theta_fail_map(:), phi_fail_map(:), ...
        p_fail_map(:), q_fail_map(:), r_fail_map(:), ...
        PLA_fail_map(:), DELTAZ_fail_map(:), DELTAW_fail_map(:), ...
        resnorm_fail_map(:), ...
        'VariableNames', { ...
        'Altitude_m','Vt_mps','TAS_mps','Feasible', ...
        'psidot_max_deg_s','phi_deg','n_max','radius_m', ...
        'alpha_deg','beta_deg','theta_deg', ...
        'p_deg_s','q_deg_s','r_deg_s', ...
        'PLA_deg','DELTAZ_deg','DELTAW_deg', ...
        'resnorm', ...
        'FailureCode','FailureReason', ...
        'psidot_fail_deg_s','alpha_fail_deg','beta_fail_deg','theta_fail_deg','phi_fail_deg', ...
        'p_fail_deg_s','q_fail_deg_s','r_fail_deg_s', ...
        'PLA_fail_deg','DELTAZ_fail_deg','DELTAW_fail_deg', ...
        'resnorm_fail'});

    writetable(T, 'tvc_turn_envelope_6dof_psidot_tas.csv');
    fprintf('Saved tvc_turn_envelope_6dof_psidot_tas.csv\n');
end

end

% =====================================================================
% SEARCH MAX PSIDOT
% =====================================================================

function [psidot_best_deg_s, y_best, info_best, ok, failInfo] = find_max_psidot_at_point( ...
    Vt, alt_m, db, opts, trim_tol, ...
    yguess_base, ...
    PLA_min, PLA_max, DELTAZ_lim, DELTAW_lim, ...
    alpha_lim, beta_lim, theta_lim_deg, phi_lim_deg, ...
    p_lim_deg_s, q_lim_deg_s, r_lim_deg_s, ...
    psidot_min_deg_s, psidot_max_deg_s, psidot_step_deg_s, psidot_tol_deg_s, ...
    vt_point_time_limit_s)

psidot_best_deg_s = NaN;
y_best = NaN(1,10);
info_best = struct('resnorm',NaN);
ok = false;

failInfo = empty_fail_info();

psidot_grid = psidot_min_deg_s:psidot_step_deg_s:psidot_max_deg_s;

last_feas_psidot = NaN;
last_feas_y = [];
last_feas_info = struct('resnorm',NaN);

first_infeas_psidot = NaN;
first_infeas_y = NaN(1,10);
first_infeas_info = struct('resnorm',NaN);
first_infeas_diag = empty_fail_info();

yguess = yguess_base;
pointClock = tic;

for k = 1:numel(psidot_grid)
    if toc(pointClock) > vt_point_time_limit_s
        failInfo.code = 15;
        failInfo.reason = "Vt-point timeout";
        failInfo.psidot_fail_deg_s = psidot_grid(k);
        failInfo.y_fail = NaN(1,10);
        failInfo.resnorm_fail = NaN;
        return
    end

    psidot_deg_s = psidot_grid(k);

    [ysol, info, feas, diag] = solve_fixed_psidot_trim( ...
        Vt, alt_m, psidot_deg_s, db, opts, trim_tol, yguess, ...
        PLA_min, PLA_max, DELTAZ_lim, DELTAW_lim, ...
        alpha_lim, beta_lim, theta_lim_deg, phi_lim_deg, ...
        p_lim_deg_s, q_lim_deg_s, r_lim_deg_s);

    if feas
        last_feas_psidot = psidot_deg_s;
        last_feas_y = ysol;
        last_feas_info = info;
        yguess = ysol;
    else
        if ~isnan(last_feas_psidot)
            first_infeas_psidot = psidot_deg_s;
            first_infeas_y = ysol;
            first_infeas_info = info;
            first_infeas_diag = diag;
            break
        else
            failInfo = diag;
            failInfo.psidot_fail_deg_s = psidot_deg_s;
            failInfo.y_fail = ysol;
            failInfo.resnorm_fail = info.resnorm;
        end
    end
end

if isnan(last_feas_psidot)
    return
end

if isnan(first_infeas_psidot)
    psidot_best_deg_s = last_feas_psidot;
    y_best = last_feas_y;
    info_best = last_feas_info;
    ok = true;

    failInfo.code = 0;
    failInfo.reason = "Reached psi_dot search ceiling";
    failInfo.psidot_fail_deg_s = psidot_max_deg_s;
    failInfo.y_fail = NaN(1,10);
    failInfo.resnorm_fail = NaN;
    return
end

psidot_lo = last_feas_psidot;
psidot_hi = first_infeas_psidot;

y_lo = last_feas_y;
info_lo = last_feas_info;

best_fail_diag = first_infeas_diag;
best_fail_y = first_infeas_y;
best_fail_info = first_infeas_info;
best_fail_psidot = first_infeas_psidot;

while (psidot_hi - psidot_lo) > psidot_tol_deg_s
    if toc(pointClock) > vt_point_time_limit_s
        failInfo = best_fail_diag;
        failInfo.code = 15;
        failInfo.reason = "Vt-point timeout";
        failInfo.psidot_fail_deg_s = best_fail_psidot;
        failInfo.y_fail = best_fail_y;
        failInfo.resnorm_fail = best_fail_info.resnorm;
        break
    end

    psidot_mid = 0.5 * (psidot_lo + psidot_hi);

    [ysol, info, feas, diag] = solve_fixed_psidot_trim( ...
        Vt, alt_m, psidot_mid, db, opts, trim_tol, y_lo, ...
        PLA_min, PLA_max, DELTAZ_lim, DELTAW_lim, ...
        alpha_lim, beta_lim, theta_lim_deg, phi_lim_deg, ...
        p_lim_deg_s, q_lim_deg_s, r_lim_deg_s);

    if feas
        psidot_lo = psidot_mid;
        y_lo = ysol;
        info_lo = info;
    else
        psidot_hi = psidot_mid;
        best_fail_diag = diag;
        best_fail_y = ysol;
        best_fail_info = info;
        best_fail_psidot = psidot_mid;
    end
end

psidot_best_deg_s = psidot_lo;
y_best = y_lo;
info_best = info_lo;
ok = true;

if failInfo.code == 0 || failInfo.reason == "Not evaluated"
    failInfo = best_fail_diag;
    failInfo.psidot_fail_deg_s = best_fail_psidot;
    failInfo.y_fail = best_fail_y;
    failInfo.resnorm_fail = best_fail_info.resnorm;
end

end

% =====================================================================
% SOLVE FIXED-PSIDOT TRIM
% =====================================================================

function [ysol, info, feasible, diag] = solve_fixed_psidot_trim( ...
    Vt, alt_m, psidot_deg_s, db, opts, trim_tol, yguess, ...
    PLA_min, PLA_max, DELTAZ_lim, DELTAW_lim, ...
    alpha_lim, beta_lim, theta_lim_deg, phi_lim_deg, ...
    p_lim_deg_s, q_lim_deg_s, r_lim_deg_s)

diag = empty_fail_info();

if isempty(yguess)
    y0 = default_turn_guess(Vt, alt_m, PLA_min, PLA_max);
else
    y0 = yguess;
end

obj = @(y) fixed_psidot_trim_residuals_safe(y, Vt, alt_m, psidot_deg_s, db);

try
    [ysol, ~, exitflag] = fsolve(obj, y0, opts);
catch
    ysol = NaN(1,10);
    info = struct('resnorm',NaN,'success',false,'exitflag',-99);
    feasible = false;
    diag.code = 2;
    diag.reason = "Solver error";
    return
end

[res, x_trim, u_trim, extra] = fixed_psidot_trim_residuals(ysol, Vt, alt_m, psidot_deg_s, db);
resnorm = norm(res);

info = struct();
info.residual = res;
info.resnorm = resnorm;
info.success = exitflag > 0 && all(isfinite(res)) && isfinite(resnorm) && resnorm < trim_tol;
info.exitflag = exitflag;
info.x_trim = x_trim;
info.u_trim = u_trim;
info.gamma_deg = rad2deg(extra.gamma_rad);
info.psi_dot_deg_s = extra.psi_dot_deg_s;
info.dx = extra.dx;

[feasible, diag] = classify_fixed_psidot_result( ...
    ysol, info, trim_tol, ...
    PLA_min, PLA_max, DELTAZ_lim, DELTAW_lim, ...
    alpha_lim, beta_lim, theta_lim_deg, phi_lim_deg, ...
    p_lim_deg_s, q_lim_deg_s, r_lim_deg_s);

end

% =====================================================================
% CLASSIFICATION
% =====================================================================

function [feasible, diag] = classify_fixed_psidot_result( ...
    ysol, info, trim_tol, ...
    PLA_min, PLA_max, DELTAZ_lim, DELTAW_lim, ...
    alpha_lim, beta_lim, theta_lim_deg, phi_lim_deg, ...
    p_lim_deg_s, q_lim_deg_s, r_lim_deg_s)

diag = empty_fail_info();
feasible = false;

if isempty(ysol) || any(~isfinite(ysol))
    diag.code = 11;
    diag.reason = "Non-finite solution";
    return
end

alpha_deg  = ysol(1);
beta_deg   = ysol(2);
theta_deg  = rad2deg(ysol(3));
phi_deg    = rad2deg(ysol(4));
p_deg_s    = rad2deg(ysol(5));
q_deg_s    = rad2deg(ysol(6));
r_deg_s    = rad2deg(ysol(7));
PLA_deg    = ysol(8);
DELTAZ_deg = ysol(9);
DELTAW_deg = ysol(10);

if info.exitflag <= 0
    diag.code = 3;
    diag.reason = "Non-converged exitflag";
    return
end

if ~(isfinite(info.resnorm) && info.resnorm < trim_tol)
    diag.code = 4;
    diag.reason = "Residual norm too large";
    return
end

if abs(alpha_deg) > alpha_lim
    diag.code = 5;
    diag.reason = "Alpha limit";
    return
end

if abs(beta_deg) > beta_lim
    diag.code = 6;
    diag.reason = "Beta limit";
    return
end

if abs(theta_deg) > theta_lim_deg
    diag.code = 7;
    diag.reason = "Theta limit";
    return
end

if abs(phi_deg) > phi_lim_deg
    diag.code = 8;
    diag.reason = "Phi limit";
    return
end

if abs(p_deg_s) > p_lim_deg_s
    diag.code = 9;
    diag.reason = "p-rate limit";
    return
end

if abs(q_deg_s) > q_lim_deg_s
    diag.code = 10;
    diag.reason = "q-rate limit";
    return
end

if abs(r_deg_s) > r_lim_deg_s
    diag.code = 11;
    diag.reason = "r-rate limit";
    return
end

if ~(PLA_deg >= PLA_min && PLA_deg <= PLA_max)
    diag.code = 12;
    diag.reason = "PLA limit";
    return
end

if abs(DELTAZ_deg) > DELTAZ_lim
    diag.code = 13;
    diag.reason = "DELTAZ limit";
    return
end

if abs(DELTAW_deg) > DELTAW_lim
    diag.code = 14;
    diag.reason = "DELTAW limit";
    return
end

feasible = true;
diag.code = 1;
diag.reason = "Feasible";

end

function diag = empty_fail_info()
diag = struct();
diag.code = 0;
diag.reason = "Not evaluated";
diag.psidot_fail_deg_s = NaN;
diag.y_fail = NaN(1,10);
diag.resnorm_fail = NaN;
end

% =====================================================================
% RESIDUALS
% =====================================================================

function R = fixed_psidot_trim_residuals_safe(y, Vt, alt_m, psidot_deg_s, db)
try
    R = fixed_psidot_trim_residuals_only(y, Vt, alt_m, psidot_deg_s, db);
    if any(~isfinite(R))
        R = 1e6 * ones(10,1);
    end
catch
    R = 1e6 * ones(10,1);
end
end

function R = fixed_psidot_trim_residuals_only(y, Vt, alt_m, psidot_deg_s, db)
[R, ~, ~, ~] = fixed_psidot_trim_residuals(y, Vt, alt_m, psidot_deg_s, db);
end

function [R, x_trim, u_trim, extra] = fixed_psidot_trim_residuals(y, Vt, alt_m, psidot_deg_s, db)
% y = [alpha_deg, beta_deg, theta_rad, phi_rad, p_rad_s, q_rad_s, r_rad_s, PLA_deg, DELTAZ_deg, DELTAW_deg]

alpha_deg  = y(1);
beta_deg   = y(2);
theta_rad  = y(3);
phi_rad    = y(4);
p_rad_s    = y(5);
q_rad_s    = y(6);
r_rad_s    = y(7);
PLA_deg    = y(8);
DELTAZ_deg = y(9);
DELTAW_deg = y(10);

if abs(alpha_deg) > 45 || ...
   abs(beta_deg) > 20 || ...
   abs(rad2deg(theta_rad)) > 60 || ...
   abs(rad2deg(phi_rad)) > 89 || ...
   abs(rad2deg(p_rad_s)) > 180 || ...
   abs(rad2deg(q_rad_s)) > 180 || ...
   abs(rad2deg(r_rad_s)) > 180 || ...
   abs(PLA_deg) > 1e4 || ...
   abs(DELTAZ_deg) > 90 || ...
   abs(DELTAW_deg) > 90

    R = 1e6 * ones(10,1);
    x_trim = nan(12,1);
    u_trim = struct('PLA',NaN,'DELTAZ',NaN,'DELTAW',NaN);
    extra = struct('gamma_rad',NaN,'psi_dot_deg_s',NaN,'dx',nan(12,1));
    return
end

alpha_rad = deg2rad(alpha_deg);
beta_rad  = deg2rad(beta_deg);

u_b = Vt * cos(alpha_rad) * cos(beta_rad);
v_b = Vt * sin(beta_rad);
w_b = Vt * sin(alpha_rad) * cos(beta_rad);

x_trim = [ ...
    u_b; ...
    v_b; ...
    w_b; ...
    p_rad_s; ...
    q_rad_s; ...
    r_rad_s; ...
    phi_rad; ...
    theta_rad; ...
    0; ...
    0; ...
    0; ...
   -alt_m];

u_trim = struct();
u_trim.PLA    = PLA_deg;
u_trim.DELTAZ = DELTAZ_deg;
u_trim.DELTAW = DELTAW_deg;

dx = ice_step_6dof(0, x_trim, u_trim, db);

phi_dot = p_rad_s + tan(theta_rad) * (q_rad_s * sin(phi_rad) + r_rad_s * cos(phi_rad));
theta_dot = q_rad_s * cos(phi_rad) - r_rad_s * sin(phi_rad);

cth = cos(theta_rad);
if abs(cth) < 1e-6
    psi_dot_model = NaN;
else
    psi_dot_model = (q_rad_s * sin(phi_rad) + r_rad_s * cos(phi_rad)) / cth;
end

gamma_rad = gamma_from_state12_local(x_trim);
psi_dot_target = deg2rad(psidot_deg_s);

R = [ ...
    dx(1); ...
    dx(2); ...
    dx(3); ...
    dx(4); ...
    dx(5); ...
    dx(6); ...
    gamma_rad; ...
    phi_dot; ...
    theta_dot; ...
    psi_dot_model - psi_dot_target ...
];

extra = struct();
extra.gamma_rad = gamma_rad;
extra.psi_dot_deg_s = rad2deg(psi_dot_model);
extra.dx = dx;

end

% =====================================================================
% INITIAL GUESS
% =====================================================================

function y0 = default_turn_guess(Vt, alt_m, PLA_min, PLA_max)
alpha0_deg  = max(1, min(10, 2 + 0.015*(200 - Vt) + 0.0001*alt_m));
beta0_deg   = 0;
theta0_rad  = deg2rad(alpha0_deg);
phi0_rad    = deg2rad(10);
p0_rad_s    = 0;
q0_rad_s    = 0;
r0_rad_s    = deg2rad(2);
PLA0_deg    = max(PLA_min + 5, min(PLA_max - 5, 50 + 0.004*alt_m + 0.04*(Vt - 100)));
DELTAZ0_deg = 0;
DELTAW0_deg = 0;

y0 = [alpha0_deg, beta0_deg, theta0_rad, phi0_rad, p0_rad_s, q0_rad_s, r0_rad_s, PLA0_deg, DELTAZ0_deg, DELTAW0_deg];
end

% =====================================================================
% GEOMETRY
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

function y = rad2deg_safe(x)
if isfinite(x)
    y = rad2deg(x);
else
    y = NaN;
end
end