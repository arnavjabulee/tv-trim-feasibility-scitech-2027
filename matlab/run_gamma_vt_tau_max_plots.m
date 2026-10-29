function results = run_gamma_vt_tau_max_plots()

db = ice_init();

% ============================================================
% USER SETTINGS
% ============================================================

alt_list = [0 10000 20000];     % [m]
Vt_list = 50:2:550;             % [m/s]
dmaxdot_list = [1 20 100 200];    % [deg/s]

delta_max_deg = 25;             % [deg]

tau_min = 0.01;                 % [s]
tau_max_cap = 10.00;            % [s]
tau_init_default = 0.50;        % [s]
tau_refine_iters = 6;
tau_expand_factor = 2.0;

tau_for_J = 1e-2;               % [s]
t_final = 20.0;                  % [s]
t_final_tau = 120.0;             % [s]

% Reject fake trims before any linearisation/simulation
trim_resnorm_max = 1e-2;

dx0 = [ ...
    0; ...
    0; ...
    0; ...
    0; ...
    deg2rad(2.0); ...
    0; ...
    0; ...
    deg2rad(10.0); ...
    0; ...
    0; ...
    0; ...
    0];

Q = diag([ ...
    1e-1 ...
    2e-1 ...
    2e-1 ...
    5e0  ...
    8e0  ...
    5e0  ...
    6e0  ...
    6e0  ...
    6e0  ...
]);

R = diag([ ...
    1e-1 ...
    5e-1 ...
]);

ode_opts = odeset('RelTol',1e-4,'AbsTol',1e-6,'MaxStep',0.05);

% ============================================================
% RUN / SAVE OPTIONS
% ============================================================

save_csv_name = 'gamma_vt_tau_max_results.csv';
save_mat_name = 'gamma_vt_tau_max_results.mat';

enable_plots = true;
enable_tau_debug_prints = false;
save_every_vpoint = true;

% ============================================================
% STORAGE
% ============================================================

results = struct();
results.alt_list = alt_list;
results.Vt_list = Vt_list;
results.dmaxdot_list = dmaxdot_list;
results.trim_resnorm_max = trim_resnorm_max;
results.tau_search = struct( ...
    'tau_min', tau_min, ...
    'tau_max_cap', tau_max_cap, ...
    'tau_init_default', tau_init_default, ...
    'tau_refine_iters', tau_refine_iters, ...
    'tau_expand_factor', tau_expand_factor);

nAlt = numel(alt_list);
nV   = numel(Vt_list);
nD   = numel(dmaxdot_list);

results.J_gamma = nan(nAlt, nD, nV);
results.tau_max = nan(nAlt, nD, nV);
results.trim_ok = false(nAlt, nV);
results.trim_resnorm = nan(nAlt, nV);

prev_tau_guess = tau_init_default * ones(nAlt, nD);

total_outer = nAlt * nV * nD;
count_outer = 0;

% ============================================================
% MAIN SWEEP
% ============================================================

for iAlt = 1:nAlt

    alt = alt_list(iAlt);

    fprintf('\n=====================================================\n');
    fprintf('Altitude = %.0f m (%d/%d)\n', alt, iAlt, nAlt);
    fprintf('=====================================================\n');

    for iV = 1:nV

        Vt = Vt_list(iV);
        fprintf('\n---- Trim / linearise at alt = %.0f m, Vt = %.1f m/s (%d/%d V-points at this altitude) ----\n', ...
            alt, Vt, iV, nV);

        % ----------------------------------------------------
        % Trim
        % ----------------------------------------------------
        try
            [x_trim12, u_trim, trim_info] = trim_6dof_level_local(Vt, alt, db);

            fprintf('  Trim candidate | alpha = %.2f deg | theta = %.2f deg | PLA = %.2f | DELTAZ = %.2f | res(trim) = %.3e | res(9) = %.3e\n', ...
                trim_info.alpha_deg, trim_info.theta_deg, trim_info.PLA_deg, ...
                trim_info.DELTAZ_deg, trim_info.resnorm_trim_eqs, trim_info.resnorm_dyn9);

        catch ME
            warning('TRIM FAILED at alt=%.0f m, Vt=%.1f m/s | %s', alt, Vt, ME.message);

            if save_every_vpoint
                save(save_mat_name, 'results');
                T_out = flatten_results_to_table(results);
                writetable(T_out, save_csv_name);
                fprintf('  Snapshot saved to %s and %s\n', save_mat_name, save_csv_name);
            end

            continue
        end

        pla_min = min(db.eng.PLA);
        pla_max = max(db.eng.PLA);

        trim_ok = ...
            all(isfinite(x_trim12)) && ...
            isfinite(u_trim.PLA) && isfinite(u_trim.DELTAZ) && isfinite(u_trim.DELTAW) && ...
            (trim_info.resnorm_trim_eqs <= trim_resnorm_max) && ...
            (u_trim.PLA >= pla_min) && (u_trim.PLA <= pla_max) && ...
            (abs(u_trim.DELTAZ) <= delta_max_deg);

        results.trim_resnorm(iAlt, iV) = trim_info.resnorm_trim_eqs;
        results.trim_ok(iAlt, iV) = trim_ok;

        if ~trim_ok
            fprintf('  Trim rejected | res(trim)=%.3e | PLA=%.2f | DELTAZ=%.2f\n', ...
                trim_info.resnorm_trim_eqs, u_trim.PLA, u_trim.DELTAZ);

            if save_every_vpoint
                save(save_mat_name, 'results');
                T_out = flatten_results_to_table(results);
                writetable(T_out, save_csv_name);
                fprintf('  Snapshot saved to %s and %s\n', save_mat_name, save_csv_name);
            end

            continue
        end

        fprintf('  Trim accepted\n');

        % ----------------------------------------------------
        % Linearise + LQR once at this trim
        % ----------------------------------------------------
        try
            [A,B] = linearize_ice_6dof(x_trim12, u_trim, db);

            Ared = A(1:9,1:9);
            Bred = B(1:9,[2 3]);

            K = lqr(Ared, Bred, Q, R);
            fprintf('  LQR OK\n');

        catch ME
            warning('LINEARISE/LQR FAILED at alt=%.0f m, Vt=%.1f m/s | %s', alt, Vt, ME.message);

            if save_every_vpoint
                save(save_mat_name, 'results');
                T_out = flatten_results_to_table(results);
                writetable(T_out, save_csv_name);
                fprintf('  Snapshot saved to %s and %s\n', save_mat_name, save_csv_name);
            end

            continue
        end

        % ----------------------------------------------------
        % For each dmaxdot:
        %   1) J(Vt) with tau ~ 0
        %   2) tau_max(Vt)
        % ----------------------------------------------------
        for iD = 1:nD

            dmaxdot = dmaxdot_list(iD);
            count_outer = count_outer + 1;

            fprintf('  [%3d/%3d] dmaxdot = %.1f deg/s\n', count_outer, total_outer, dmaxdot);

            % ================================================
            % 1) J(Vt) at tau ~ 0
            % ================================================
            try
                perfJ = default_tvc_perf();
                perfJ.delta_max_deg = delta_max_deg;
                perfJ.delta_dot_max_deg_s = dmaxdot;
                perfJ.tau_tvc_s = tau_for_J;
                perfJ.tau_thrust_s = 0.01;
                perfJ.PLA_dot_max_per_s = inf;

                x0_12 = x_trim12 + dx0;
                x0_15 = [x0_12; u_trim.DELTAZ; u_trim.DELTAW; u_trim.PLA];
                x0_aug = [x0_15; 0];

                rhs = @(t,x) rhs_trim_regulation_local(t, x, x_trim12, u_trim, K, db, perfJ);

                fprintf('    starting J sim: tau=%.4f, Vt=%.1f, alt=%.0f\n', perfJ.tau_tvc_s, Vt, alt);
                tic
                [t_hist, x_aug_hist] = ode45(rhs, [0 t_final], x0_aug, ode_opts);
                fprintf('    finished J sim in %.2f s\n', toc);

                x15_hist = x_aug_hist(:,1:15);
                x12_hist = x15_hist(:,1:12);
                [~,~,~,gamma_deg] = derive_airdata_gamma_from_states(x12_hist);

                J = trapz(t_hist, abs(gamma_deg));
                results.J_gamma(iAlt, iD, iV) = J;

                fprintf('    J(tau≈0) = %.3f deg*s\n', J);

            catch ME
                warning('J SIM FAILED at alt=%.0f m, Vt=%.1f m/s, dmaxdot=%.1f | %s', ...
                    alt, Vt, dmaxdot, ME.message);
            end

            % ================================================
            % 2) tau_max(Vt)
            % ================================================
            tau_guess = prev_tau_guess(iAlt, iD);
            tau_guess = min(max(tau_guess, tau_min), tau_max_cap);

            fprintf('    tau search starts from guess = %.3f s\n', tau_guess);

            [stable_guess, ~] = run_tau_stability_test_local( ...
                tau_guess, x_trim12, u_trim, K, db, delta_max_deg, dmaxdot, ...
                dx0, t_final_tau, ode_opts, enable_tau_debug_prints);

            if stable_guess
                tau_lo = tau_guess;
                tau_hi = NaN;

                tau_try = tau_guess;
                while tau_try < tau_max_cap
                    tau_try = min(tau_try * tau_expand_factor, tau_max_cap);

                    [stable_try, ~] = run_tau_stability_test_local( ...
                        tau_try, x_trim12, u_trim, K, db, delta_max_deg, dmaxdot, ...
                        dx0, t_final_tau, ode_opts, enable_tau_debug_prints);

                    if stable_try
                        tau_lo = tau_try;
                        if abs(tau_try - tau_max_cap) < 1e-12
                            break
                        end
                    else
                        tau_hi = tau_try;
                        break
                    end
                end

                if isnan(tau_hi)
                    tau_max_this = tau_lo;
                else
                    for ib = 1:tau_refine_iters
                        tau_mid = 0.5 * (tau_lo + tau_hi);

                        [stable_mid, ~] = run_tau_stability_test_local( ...
                            tau_mid, x_trim12, u_trim, K, db, delta_max_deg, dmaxdot, ...
                            dx0, t_final_tau, ode_opts, enable_tau_debug_prints);

                        if stable_mid
                            tau_lo = tau_mid;
                        else
                            tau_hi = tau_mid;
                        end
                    end
                    tau_max_this = tau_lo;
                end

            else
                tau_hi = tau_guess;
                tau_lo = NaN;

                tau_try = tau_guess;
                while tau_try > tau_min
                    tau_try = max(tau_try / tau_expand_factor, tau_min);

                    [stable_try, ~] = run_tau_stability_test_local( ...
                        tau_try, x_trim12, u_trim, K, db, delta_max_deg, dmaxdot, ...
                        dx0, t_final_tau, ode_opts, enable_tau_debug_prints);

                    if stable_try
                        tau_lo = tau_try;
                        break
                    elseif abs(tau_try - tau_min) < 1e-12
                        break
                    end
                end

                if isnan(tau_lo)
                    tau_max_this = NaN;
                else
                    for ib = 1:tau_refine_iters
                        tau_mid = 0.5 * (tau_lo + tau_hi);

                        [stable_mid, ~] = run_tau_stability_test_local( ...
                            tau_mid, x_trim12, u_trim, K, db, delta_max_deg, dmaxdot, ...
                            dx0, t_final_tau, ode_opts, enable_tau_debug_prints);

                        if stable_mid
                            tau_lo = tau_mid;
                        else
                            tau_hi = tau_mid;
                        end
                    end
                    tau_max_this = tau_lo;
                end
            end

            results.tau_max(iAlt, iD, iV) = tau_max_this;

            if ~isnan(tau_max_this)
                prev_tau_guess(iAlt, iD) = tau_max_this;
            end

            fprintf('    tau_max = %.3f s\n', tau_max_this);
        end

        if save_every_vpoint
            save(save_mat_name, 'results');
            T_out = flatten_results_to_table(results);
            writetable(T_out, save_csv_name);
            fprintf('  Snapshot saved to %s and %s\n', save_mat_name, save_csv_name);
        end
    end
end

save(save_mat_name, 'results');
T_out = flatten_results_to_table(results);
writetable(T_out, save_csv_name);

if enable_plots
    plot_results_local(results);
end

end

% ============================================================
% PLOTTING
% ============================================================

function plot_results_local(results)

nr = numel(results.dmaxdot_list);

labels = strings(nr,1);
for i = 1:nr
    labels(i) = sprintf('\\delta-dot_{max} = %.1f deg/s', results.dmaxdot_list(i));
end

C = lines(nr);

for iAlt = 1:numel(results.alt_list)

    alt = results.alt_list(iAlt);

    figure('Name', sprintf('Gamma tracking summary | h = %.0f m', alt));
    tiledlayout(1,2);

    nexttile;
    hold on
    for iD = 1:nr
        y = squeeze(results.J_gamma(iAlt, iD, :));
        plot(results.Vt_list, y, '-o', ...
            'LineWidth', 1.8, ...
            'MarkerSize', 5, ...
            'Color', C(iD,:));
    end
    grid on
    xlabel('V_t (m/s)')
    ylabel('J = IAE_\gamma (deg*s)')
    title(sprintf('Performance at \\tau \\approx 0 | h = %.0f m', alt))

    nexttile;
    hold on
    for iD = 1:nr
        y = squeeze(results.tau_max(iAlt, iD, :));
        plot(results.Vt_list, y, '-o', ...
            'LineWidth', 1.8, ...
            'MarkerSize', 5, ...
            'Color', C(iD,:));
    end
    grid on
    xlabel('V_t (m/s)')
    ylabel('\tau_{max} (s)')
    title(sprintf('Maximum stable actuator lag | h = %.0f m', alt))

    legend(labels,'Location','northoutside','Orientation','horizontal')
end

end

% ============================================================
% TRIM
% ============================================================

function [x_trim, u_trim, info] = trim_6dof_level_local(Vt_target, alt_m, db)

zD = -alt_m;

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-5, ...
    'StepTolerance',1e-5, ...
    'MaxFunctionEvaluations',800, ...
    'MaxIterations',100);

x0 = [5; deg2rad(5); 87; 0];
sol = fsolve(@(y) trim_residuals_local(y, Vt_target, zD, db), x0, opts);

alpha_deg  = sol(1);
theta_rad  = sol(2);
PLA_deg    = sol(3);
DELTAZ_deg = sol(4);

alpha_rad = deg2rad(alpha_deg);

u = Vt_target * cos(alpha_rad);
v = 0;
w = Vt_target * sin(alpha_rad);

x_trim = [u; v; w; 0; 0; 0; 0; theta_rad; 0; 0; 0; zD];

u_trim = struct();
u_trim.PLA = PLA_deg;
u_trim.DELTAZ = DELTAZ_deg;
u_trim.DELTAW = 0;

dx = ice_step_6dof(0, x_trim, u_trim, db);

info = struct();
info.alpha_deg = alpha_deg;
info.theta_deg = rad2deg(theta_rad);
info.PLA_deg = PLA_deg;
info.DELTAZ_deg = DELTAZ_deg;
info.dx_trim = dx;
info.residual_trim_eqs = [dx(1); dx(3); dx(5); dx(12)];
info.resnorm_trim_eqs = norm(info.residual_trim_eqs);
info.resnorm_dyn9 = norm(dx(1:9));

end

function R = trim_residuals_local(y, Vt_target, zD, db)

alpha_deg  = y(1);
theta_rad  = y(2);
PLA_deg    = y(3);
DELTAZ_deg = y(4);

alpha_rad = deg2rad(alpha_deg);

u = Vt_target * cos(alpha_rad);
v = 0;
w = Vt_target * sin(alpha_rad);

x = [u; v; w; 0; 0; 0; 0; theta_rad; 0; 0; 0; zD];

u_in = struct();
u_in.PLA = PLA_deg;
u_in.DELTAZ = DELTAZ_deg;
u_in.DELTAW = 0;

dx = ice_step_6dof(0, x, u_in, db);

R = [dx(1); dx(3); dx(5); dx(12)];

end

% ============================================================
% CLOSED-LOOP RHS
% 16 states = 15 plant + I_dummy
% ============================================================

function dx_aug = rhs_trim_regulation_local(~, x_aug, x_trim12, u_trim, K, db, perf)

x15 = x_aug(1:15);
x12 = x15(1:12);

xd = x12(1:9) - x_trim12(1:9);
du = -K * xd;

u = u_trim;
u.DELTAZ = saturate_deg_local(u_trim.DELTAZ + du(1), perf.delta_max_deg);
u.DELTAW = saturate_deg_local(u_trim.DELTAW + du(2), perf.delta_max_deg);
u.PLA = u_trim.PLA;

dx15 = ice_step_6dof(0, x15, u, db, perf);
I_dummy_dot = 0;

dx_aug = [dx15; I_dummy_dot];

end

% ============================================================
% STABILITY TEST
% ============================================================

function [is_stable, info] = run_tau_stability_test_local( ...
    tau, x_trim12, u_trim, K, db, delta_max_deg, dmaxdot, dx0, t_final_tau, ode_opts, enable_tau_debug_prints)

perfTau = default_tvc_perf();
perfTau.delta_max_deg = delta_max_deg;
perfTau.delta_dot_max_deg_s = dmaxdot;
perfTau.tau_tvc_s = max(tau, 1e-2);
perfTau.tau_thrust_s = 0.01;
perfTau.PLA_dot_max_per_s = inf;

x0_12 = x_trim12 + dx0;
x0_15 = [x0_12; u_trim.DELTAZ; u_trim.DELTAW; u_trim.PLA];
x0_aug = [x0_15; 0];

rhsTau = @(t,x) rhs_trim_regulation_local(t, x, x_trim12, u_trim, K, db, perfTau);

fprintf('    starting tau test: tau=%.3f s\n', tau);
tic
try
    [t_tau, x_tau] = ode45(rhsTau, [0 t_final_tau], x0_aug, ode_opts);
    ode_fail = false;
    ode_fail_msg = '';
catch ME
    t_tau = [];
    x_tau = [];
    ode_fail = true;
    ode_fail_msg = ME.message;
end
fprintf('    finished tau test in %.2f s\n', toc);

[sim_ok, fail_reason] = sim_health_check_local(t_tau, x_tau, t_final_tau, ode_fail, ode_fail_msg);

if isempty(t_tau) || isempty(x_tau)
    is_stable = false;
    info = struct();
    info.t = t_tau;
    info.x = x_tau;
    info.gamma = [];
    info.max_gamma = NaN;
    info.final_gamma = NaN;
    info.conv_info = [];
    info.sim_ok = sim_ok;
    info.fail_reason = fail_reason;
    fprintf('    tau = %.3f s | stable = 0 | sim failed: %s\n', tau, fail_reason);
    return
end

x15_tau = x_tau(:,1:15);

% ------------------------------------------------------------
% DEBUG: check actuator rate usage
% ------------------------------------------------------------
DELTAZ_act = x15_tau(:,13);
DELTAW_act = x15_tau(:,14);

dt = diff(t_tau);

if enable_tau_debug_prints && ~isempty(dt)
    delz_rate_est = diff(DELTAZ_act) ./ dt;
    delw_rate_est = diff(DELTAW_act) ./ dt;

    fprintf('    max |DELTAZ_dot| = %.3f deg/s | limit = %.3f\n', ...
        max(abs(delz_rate_est)), dmaxdot);

    fprintf('    max |DELTAW_dot| = %.3f deg/s | limit = %.3f\n', ...
        max(abs(delw_rate_est)), dmaxdot);
end

x12_tau = x15_tau(:,1:12);

[~,~,~,gamma_tau] = derive_airdata_gamma_from_states(x12_tau);
psi_tau_deg = rad2deg(x12_tau(:,9));

max_gamma = max(abs(gamma_tau));
final_gamma = abs(gamma_tau(end));

if sim_ok
    [is_converging, conv_info] = trim_regulation_convergence_test_local(t_tau, gamma_tau, 0.0, psi_tau_deg);
else
    is_converging = false;
    conv_info = default_conv_info_failed();
end

is_stable = sim_ok && max_gamma < 45 && is_converging;

info = struct();
info.t = t_tau;
info.x = x_tau;
info.gamma = gamma_tau;
info.max_gamma = max_gamma;
info.final_gamma = final_gamma;
info.conv_info = conv_info;
info.sim_ok = sim_ok;
info.fail_reason = fail_reason;

fprintf('    tau = %.3f s | stable = %d | method = %s | max|gamma| = %.3f | final|gamma| = %.3f\n', ...
    tau, is_stable, conv_info.method, max_gamma, final_gamma);

end

% ============================================================
% CONVERGENCE TEST
% ============================================================

function [is_converging, info] = trim_regulation_convergence_test_local(t, y, y_cmd, psi_deg)

t = t(:);
y = y(:);
psi_deg = psi_deg(:);

info = struct();
info.extrema_idx = [];
info.amp_seq = [];
info.method = "invalid";
info.tail_slope = NaN;
info.tail_mean_1 = NaN;
info.tail_mean_2 = NaN;
info.final_abs_error = NaN;
info.start_idx = NaN;
info.psi_fail = false;
info.psi_abs_max = NaN;
info.psi_tail_mean_1 = NaN;
info.psi_tail_mean_2 = NaN;
info.psi_tail_slope = NaN;

if numel(t) ~= numel(y) || numel(y) < 6 || any(~isfinite(t)) || any(~isfinite(y)) || ...
        numel(psi_deg) ~= numel(t) || any(~isfinite(psi_deg))
    is_converging = false;
    return
end

e = y - y_cmd;
ae = abs(e);
apsi = abs(psi_deg);

t_ignore = 2.0;
min_tail_pts = 12;
tail_decay_ratio = 0.85;
final_error_tol = 2.0;
amp_shrink_ratio = 0.85;

psi_abs_limit = 20.0;
psi_tail_growth_ratio = 1.10;

i_start = find(t >= t_ignore, 1, 'first');
if isempty(i_start)
    i_start = 1;
end
info.start_idx = i_start;

t_use = t(i_start:end);
e_use = e(i_start:end);
ae_use = ae(i_start:end);
apsi_use = apsi(i_start:end);

if numel(t_use) < min_tail_pts
    is_converging = false;
    info.method = "too_short";
    return
end

info.psi_abs_max = max(apsi_use);

nPsi = numel(apsi_use);
halfPsi = floor(nPsi/2);

if halfPsi >= 3 && (nPsi - halfPsi) >= 3
    psi_tail_1 = apsi_use(1:halfPsi);
    psi_tail_2 = apsi_use(halfPsi+1:end);
    info.psi_tail_mean_1 = mean(psi_tail_1);
    info.psi_tail_mean_2 = mean(psi_tail_2);
    ppsi = polyfit(t_use, apsi_use, 1);
    info.psi_tail_slope = ppsi(1);
end

psi_hard_fail = info.psi_abs_max > psi_abs_limit;
psi_tail_fail = ...
    isfinite(info.psi_tail_mean_1) && ...
    isfinite(info.psi_tail_mean_2) && ...
    (info.psi_tail_mean_2 > psi_tail_growth_ratio * info.psi_tail_mean_1) && ...
    isfinite(info.psi_tail_slope) && ...
    (info.psi_tail_slope > 0);

if psi_hard_fail || psi_tail_fail
    is_converging = false;
    info.psi_fail = true;
    info.method = "psi_divergence";
    info.final_abs_error = ae(end);
    return
end

de = diff(e_use) ./ max(diff(t_use), 1e-12);
s = sign(de);

for k = 2:numel(s)
    if s(k) == 0
        s(k) = s(k-1);
    end
end
for k = numel(s)-1:-1:1
    if s(k) == 0
        s(k) = s(k+1);
    end
end

extrema_idx_local = [];
for k = 1:numel(s)-1
    if s(k) ~= 0 && s(k+1) ~= 0 && s(k+1) ~= s(k)
        extrema_idx_local(end+1,1) = k+1; %#ok<AGROW>
    end
end

if ~isempty(extrema_idx_local)
    keep = true(size(extrema_idx_local));
    for k = 2:numel(extrema_idx_local)
        if extrema_idx_local(k) - extrema_idx_local(k-1) <= 1
            keep(k) = false;
        end
    end
    extrema_idx_local = extrema_idx_local(keep);
end

amp_seq = abs(e_use(extrema_idx_local));
extrema_idx = extrema_idx_local + i_start - 1;

info.extrema_idx = extrema_idx;
info.amp_seq = amp_seq;
info.final_abs_error = ae(end);

if numel(amp_seq) >= 2
    overall_shrink = amp_seq(end) < amp_shrink_ratio * amp_seq(1);
    kk = (1:numel(amp_seq)).';
    p_ext = polyfit(kk, amp_seq, 1);
    slope_negative = p_ext(1) < 0;
    ratio = amp_seq(2:end) ./ max(amp_seq(1:end-1), 1e-12);
    n_bad_jumps = sum(ratio > 1.10);
    acceptable_local_behaviour = n_bad_jumps <= max(1, floor(0.25 * numel(ratio)));

    is_converging = overall_shrink && slope_negative && acceptable_local_behaviour;
    info.method = "extrema";
    return
end

nTail = numel(ae_use);
halfIdx = floor(nTail/2);

if halfIdx < 3 || (nTail - halfIdx) < 3
    is_converging = false;
    info.method = "tail_too_short";
    return
end

tail_1 = ae_use(1:halfIdx);
tail_2 = ae_use(halfIdx+1:end);

tail_mean_1 = mean(tail_1);
tail_mean_2 = mean(tail_2);

p_tail = polyfit(t_use, ae_use, 1);
tail_slope = p_tail(1);

info.tail_mean_1 = tail_mean_1;
info.tail_mean_2 = tail_mean_2;
info.tail_slope = tail_slope;
info.method = "windowed_tail";

is_converging = ...
    isfinite(tail_mean_1) && isfinite(tail_mean_2) && isfinite(tail_slope) && ...
    (tail_mean_2 < tail_decay_ratio * tail_mean_1) && ...
    (tail_slope < 0) && ...
    (ae(end) < final_error_tol);

end

function info = default_conv_info_failed()
info = struct();
info.method = "sim_failed";
info.final_abs_error = NaN;
info.amp_seq = [];
info.extrema_idx = [];
info.tail_mean_1 = NaN;
info.tail_mean_2 = NaN;
info.tail_slope = NaN;
info.psi_fail = false;
info.psi_abs_max = NaN;
info.psi_tail_mean_1 = NaN;
info.psi_tail_mean_2 = NaN;
info.psi_tail_slope = NaN;
end

% ============================================================
% AIRDATA / GAMMA
% ============================================================

function [Vt, alpha_deg, beta_deg, gamma_deg] = derive_airdata_gamma_from_states(xhist)

u = xhist(:,1);
v = xhist(:,2);
w = xhist(:,3);

phi   = xhist(:,7);
theta = xhist(:,8);
psi   = xhist(:,9);

Vt = sqrt(u.^2 + v.^2 + w.^2);
alpha_deg = rad2deg(atan2(w, u));
beta_deg = rad2deg(asin(max(min(v ./ max(Vt,1e-9), 1), -1)));

n = size(xhist,1);
gamma_deg = zeros(n,1);

for i = 1:n
    cphi = cos(phi(i));   sphi = sin(phi(i));
    cth  = cos(theta(i)); sth  = sin(theta(i));
    cpsi = cos(psi(i));   spsi = sin(psi(i));

    R_b2n = [ cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
              cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
                 -sth,                sphi*cth,                 cphi*cth      ];

    Vn = R_b2n * [u(i); v(i); w(i)];
    gamma_deg(i) = rad2deg(atan2(-Vn(3), sqrt(Vn(1)^2 + Vn(2)^2)));
end

end

% ============================================================
% SIM HEALTH
% ============================================================

function [sim_ok, fail_reason] = sim_health_check_local(t, x_aug, tfinal, ode_fail, ode_fail_msg)

sim_ok = true;
fail_reason = 'ok';

if nargin >= 4 && ode_fail
    sim_ok = false;
    fail_reason = ['ode45_failed: ' ode_fail_msg];
    return
end

if isempty(t) || isempty(x_aug)
    sim_ok = false;
    fail_reason = 'empty_output';
    return
end

if any(~isfinite(t)) || any(~isfinite(x_aug(:)))
    sim_ok = false;
    fail_reason = 'nonfinite';
    return
end

if t(end) < tfinal - 1e-6
    sim_ok = false;
    fail_reason = 'terminated_early';
    return
end

end

% ============================================================
% CSV FLATTENING
% ============================================================

function T = flatten_results_to_table(results)

alt_list = results.alt_list;
Vt_list = results.Vt_list;
dmaxdot_list = results.dmaxdot_list;

nAlt = numel(alt_list);
nV = numel(Vt_list);
nD = numel(dmaxdot_list);

Alt_m = [];
Vt_mps = [];
DeltaDotMax_degps = [];
TrimOK = [];
TrimResnorm = [];
J_gamma_deg_s = [];
TauMax_s = [];

for iAlt = 1:nAlt
    for iV = 1:nV
        for iD = 1:nD
            Alt_m(end+1,1) = alt_list(iAlt); %#ok<AGROW>
            Vt_mps(end+1,1) = Vt_list(iV); %#ok<AGROW>
            DeltaDotMax_degps(end+1,1) = dmaxdot_list(iD); %#ok<AGROW>

            TrimOK(end+1,1) = results.trim_ok(iAlt, iV); %#ok<AGROW>
            TrimResnorm(end+1,1) = results.trim_resnorm(iAlt, iV); %#ok<AGROW>

            J_gamma_deg_s(end+1,1) = results.J_gamma(iAlt, iD, iV); %#ok<AGROW>
            TauMax_s(end+1,1) = results.tau_max(iAlt, iD, iV); %#ok<AGROW>
        end
    end
end

T = table( ...
    Alt_m, ...
    Vt_mps, ...
    DeltaDotMax_degps, ...
    TrimOK, ...
    TrimResnorm, ...
    J_gamma_deg_s, ...
    TauMax_s);

end

% ============================================================
% SATURATION HELPERS
% ============================================================

function y = saturate_deg_local(x, lim)
y = min(max(x, -lim), lim);
end