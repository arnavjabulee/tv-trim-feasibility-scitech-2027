function sims = run_lqr_tvc_cases()
% RUN_LQR_TVC_CASES
%
% TVC-only trim-regulation study:
%   - exact steady-level 6DOF trim
%   - linearise at trim
%   - LQR on reduced 10-state plant using TVC only [DELTAZ DELTAW]
%   - reduced controller states:
%         [u v w p q r phi theta psi zD]
%   - PLA held fixed at PLA_trim
%   - nonlinear simulation uses 15-state plant:
%         x15 = [12 aircraft states; DELTAZ_act; DELTAW_act; PLA_act]
%
% Requires:
%   ice_init.m
%   ice_step_6dof.m
%   linearize_ice_6dof.m
%   default_tvc_perf.m
%   Control System Toolbox

db = ice_init();

% ============================================================
% EDIT CASES HERE
% dx0 is a 12x1 perturbation on:
% [u v w p q r phi theta psi xN yE zD]
% ============================================================
cases = [
    make_case("main envelope", 250, 6000, ...
       [0;0;0; 0;deg2rad(2.0);0; 0;deg2rad(10.0);0; 0;0;0], 120)

    % make_case("island", 480, 18000, ...
    %    [0;0;0; 0;deg2rad(2.0);0; 0;deg2rad(10.0);0; 0;0;0], 120)

    % make_case("unstable 1", 450, 16000, ...
    %    [0;0;0; 0;deg2rad(2.0);0; 0;deg2rad(10.0);0; 0;0;0], 120)

    % make_case("unstable 2", 200, 14000, ...
    %    [0;0;0; 0;deg2rad(2.0);0; 0;deg2rad(0.0);0; 0;0;0], 120)

    % make_case("unstable 3", 450, 6000, ...
    %    [0;0;0; 0;deg2rad(2.0);0; 0;deg2rad(0.0);0; 0;0;0], 120)

    % make_case("Low-speed check", 76, 0, ...
    %    [0;0;0; 0;deg2rad(2.0);0; 0;deg2rad(10.0);0; 0;0;0], 20)
];

% ============================================================
% TVC / engine dynamics performance model
% ============================================================
perf = default_tvc_perf();
perf.delta_max_deg = 25;
perf.delta_dot_max_deg_s = 50;
perf.tau_tvc_s = 0.1;
perf.tau_thrust_s = 0.42;
perf.PLA_dot_max_per_s = inf;

% Controller design state indices:
% [u v w p q r phi theta psi zD]
ix_lqr = [1 2 3 4 5 6 7 8 9 12];

% LQR weights for reduced 10-state controller design
Q_lqr = diag([ ...
    5e0 ... % u      lower than before
    1e-1 ... % v
    1e-1 ... % w
    2e0  ... % p      lower rate penalty
    3e0  ... % q      lower q penalty
    2e0  ... % r
    3e0  ... % phi
    3e0  ... % theta
    5e-1 ... % psi
    5e-1 ... % zD     much lower altitude/position penalty
]);

R = diag([ ...
    2e0 ... % DELTAZ, was 1e-1
    5e0 ... % DELTAW, was 5e-1
]);


ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.02);

sims = struct([]);

for k = 1:numel(cases)
    c = cases(k);

    fprintf('\n================ %s ================\n', c.name);
    fprintf('Target: Vt = %.1f m/s, alt = %.1f m\n', c.Vt, c.alt_m);

    % --------------------------------------------------------
    % Trim exact 6-DOF point
    % --------------------------------------------------------
    [x_trim12, u_trim, info] = trim_6dof_level_local(c.Vt, c.alt_m, db);

    pla_margin_deg = 1.0;
    gamma_tol_deg = 0.1;
    pla_max_deg = 127.0;

    trim_ok = ...
        isfinite(info.resnorm) && ...
        info.resnorm <= 1e-5 && ...
        info.exitflag > 0 && ...
        abs(info.gamma_deg) <= gamma_tol_deg && ...
        (info.PLA_deg < pla_max_deg - pla_margin_deg);

    fprintf('Trim:\n');
    fprintf('  alpha  = %.3f deg\n', info.alpha_deg);
    fprintf('  theta  = %.3f deg\n', info.theta_deg);
    fprintf('  gamma  = %.6f deg\n', info.gamma_deg);
    fprintf('  PLA    = %.3f deg\n', info.PLA_deg);
    fprintf('  DELTAZ = %.3f deg\n', info.DELTAZ_deg);
    fprintf('  resnorm= %.3e\n', info.resnorm);

    if ~trim_ok
        warning('Skipping case "%s": trim failed. exitflag = %g, resnorm = %.3e, PLA = %.3f deg, gamma = %.3f deg.', ...
            c.name, info.exitflag, info.resnorm, info.PLA_deg, info.gamma_deg);

        sims(k) = build_failed_case_struct( ...
            c, x_trim12, u_trim, info, [], [], [], perf, ...
            'trim_failed', false);

        save_case_csv_local(sims(k));
        continue
    end

    % --------------------------------------------------------
    % Linearise
    % --------------------------------------------------------
    [A,B] = linearize_ice_6dof(x_trim12, u_trim, db);

    % Reduced controller-design plant:
    % [u v w p q r phi theta psi zD]
    A10 = A(ix_lqr, ix_lqr);

    % Assumed input ordering in linearize_ice_6dof:
    % [PLA DELTAZ DELTAW]
    B10 = B(ix_lqr, [2 3]);

    fprintf('rank(ctrb(A10,B10)) = %d / %d\n', ...
        rank(ctrb(A10,B10)), size(A10,1));

    fprintf('Open-loop eigenvalues, A10:\n');
    disp(eig(A10));

    % --------------------------------------------------------
    % LQR, no integrator
    % --------------------------------------------------------
    Kx = lqr(A10, B10, Q_lqr, R);

    fprintf('Closed-loop eigenvalues, A10 - B10*Kx:\n');
    disp(eig(A10 - B10*Kx));

    % --------------------------------------------------------
    % Simulate
    % state = x15 = [12 aircraft; DELTAZ_act; DELTAW_act; PLA_act]
    % --------------------------------------------------------
    x0_12 = x_trim12 + c.dx0;
    x0_15 = [x0_12; u_trim.DELTAZ; u_trim.DELTAW; u_trim.PLA];

    rhs = @(t,x) rhs_trim_regulation_local(t, x, x_trim12, u_trim, Kx, db, perf, ix_lqr);

    fprintf('Starting nonlinear sim...\n');
    tic
    try
        [t_cl, x15_hist] = ode45(rhs, [0 c.tfinal], x0_15, ode_opts);
        ode_fail = false;
        ode_fail_msg = '';
    catch ME
        t_cl = [];
        x15_hist = [];
        ode_fail = true;
        ode_fail_msg = ME.message;
    end
    sim_time = toc;

    [sim_ok, fail_reason] = sim_health_check_local(t_cl, x15_hist, c.tfinal, ode_fail, ode_fail_msg);

    if isempty(t_cl)
        t_end = NaN;
    else
        t_end = t_cl(end);
    end

    fprintf('Simulation health: ok = %d | reason = %s | tfinal reached = %.3f / %.3f\n', ...
        sim_ok, fail_reason, t_end, c.tfinal);
    fprintf('Finished in %.2f s\n', sim_time);

    if isempty(t_cl) || isempty(x15_hist)
        sims(k) = build_failed_case_struct(c, x_trim12, u_trim, info, A, B, Kx, perf, fail_reason, trim_ok);
        save_case_csv_local(sims(k));
        continue
    end

    x_cl = x15_hist(:,1:12);

    DELTAZ_act = x15_hist(:,13);
    DELTAW_act = x15_hist(:,14);
    PLA_act    = x15_hist(:,15);

    [Vt_cl, alpha_cl, beta_cl, gamma_cl] = derive_airdata_gamma_from_states(x_cl);

    [delz_cmd, delw_cmd, pla_cmd] = recover_trim_regulation_signals( ...
        x15_hist, x_trim12, u_trim, Kx, perf, ix_lqr);

    psi_cl_deg = rad2deg(x_cl(:,9));
    alt_cl = -x_cl(:,12);
    alt_trim = -x_trim12(12);

    if sim_ok
        [is_converging, conv_info] = trim_regulation_convergence_test_local(t_cl, alt_cl, alt_trim, psi_cl_deg);
        amp_seq = conv_info.amp_seq;
        extrema_idx = conv_info.extrema_idx;
    else
        is_converging = false;
        conv_info = struct();
        conv_info.method = "sim_failed";
        conv_info.final_abs_error = NaN;
        conv_info.amp_seq = [];
        conv_info.extrema_idx = [];
        conv_info.tail_mean_1 = NaN;
        conv_info.tail_mean_2 = NaN;
        conv_info.tail_slope = NaN;
        conv_info.psi_fail = false;
        conv_info.psi_abs_max = NaN;
        conv_info.psi_tail_mean_1 = NaN;
        conv_info.psi_tail_mean_2 = NaN;
        conv_info.psi_tail_slope = NaN;
        amp_seq = [];
        extrema_idx = [];
    end

    fprintf('Altitude convergence test:\n');
    fprintf('  converging = %d | method = %s | final |h-h_trim| = %.3f m\n', ...
        is_converging, conv_info.method, conv_info.final_abs_error);

    if conv_info.psi_fail
        fprintf('  PSI FAIL | max|psi| = %.3f | psi tail mean 1 = %.3f | psi tail mean 2 = %.3f | psi tail slope = %.5f\n', ...
            conv_info.psi_abs_max, conv_info.psi_tail_mean_1, conv_info.psi_tail_mean_2, conv_info.psi_tail_slope);
    elseif strcmp(conv_info.method, "extrema")
        fprintf('  first amp = %.3f | last amp = %.3f | nExt = %d\n', ...
            conv_info.amp_seq(1), conv_info.amp_seq(end), numel(conv_info.extrema_idx));
    elseif strcmp(conv_info.method, "windowed_tail")
        fprintf('  alt tail mean 1 = %.3f | alt tail mean 2 = %.3f | alt tail slope = %.5f\n', ...
            conv_info.tail_mean_1, conv_info.tail_mean_2, conv_info.tail_slope);
    end

    sims(k).name = c.name;
    sims(k).case = c;
    sims(k).x_trim12 = x_trim12;
    sims(k).u_trim = u_trim;
    sims(k).trim_info = info;
    sims(k).A = A;
    sims(k).B = B;
    sims(k).Kx = Kx;
    sims(k).K_lqr = Kx;
    sims(k).perf = perf;
    sims(k).conv_info = conv_info;

    sims(k).trim_ok = trim_ok;
    sims(k).sim_ok = sim_ok;
    sims(k).fail_reason = fail_reason;

    sims(k).t_cl = t_cl;
    sims(k).x15 = x15_hist;
    sims(k).x_cl = x_cl;

    sims(k).Vt_cl = Vt_cl;
    sims(k).alpha_cl = alpha_cl;
    sims(k).beta_cl = beta_cl;
    sims(k).gamma_cl = gamma_cl;
    sims(k).alt_cl = alt_cl;
    sims(k).alt_trim = alt_trim;

    sims(k).delz_cmd = delz_cmd;
    sims(k).delw_cmd = delw_cmd;
    sims(k).pla_cmd = pla_cmd;

    sims(k).DELTAZ_act = DELTAZ_act;
    sims(k).DELTAW_act = DELTAW_act;
    sims(k).PLA_act = PLA_act;

    sims(k).is_converging = is_converging;
    sims(k).amp_seq = amp_seq;
    sims(k).extrema_idx = extrema_idx;

    plot_case_results_local(sims(k));
    save_case_csv_local(sims(k));
end

end

% ============================================================
% FAILED CASE BUILDER
% ============================================================

function s = build_failed_case_struct(c, x_trim12, u_trim, info, A, B, Kx, perf, fail_reason, trim_ok)

s.name = c.name;
s.case = c;
s.x_trim12 = x_trim12;
s.u_trim = u_trim;
s.trim_info = info;
s.A = A;
s.B = B;
s.Kx = Kx;
s.K_lqr = Kx;
s.perf = perf;

s.trim_ok = trim_ok;
s.sim_ok = false;
s.fail_reason = fail_reason;

s.t_cl = [];
s.x15 = [];
s.x_cl = [];

s.Vt_cl = [];
s.alpha_cl = [];
s.beta_cl = [];
s.gamma_cl = [];
s.alt_cl = [];
s.alt_trim = [];

s.delz_cmd = [];
s.delw_cmd = [];
s.pla_cmd = [];

s.DELTAZ_act = [];
s.DELTAW_act = [];
s.PLA_act = [];

s.is_converging = false;
s.amp_seq = [];
s.extrema_idx = [];

end

% ============================================================
% CASE HELPER
% ============================================================

function c = make_case(name, Vt, alt_m, dx0, tfinal)
c = struct();
c.name = name;
c.Vt = Vt;
c.alt_m = alt_m;
c.dx0 = dx0(:);
c.tfinal = tfinal;
end

% ============================================================
% TRIM
% ============================================================

function [x_trim, u_trim, info] = trim_6dof_level_local(Vt_target, alt_m, db)

zD = -alt_m;

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-8, ...
    'MaxFunctionEvaluations',2000, ...
    'MaxIterations',250);

% unknowns = [alpha_deg, theta_rad, PLA_deg, DELTAZ_deg]
x0 = [5; deg2rad(5); 87; 0];

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

if abs(alpha_deg) > 40 || ...
   abs(rad2deg(theta_rad)) > 50 || ...
   PLA_deg < -20 || PLA_deg > 127 || ...
   abs(DELTAZ_deg) > 50

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
    gamma_actual_rad ...
];

if any(~isfinite(R))
    R = 1e6 * ones(4,1);
end

end

% ============================================================
% CLOSED-LOOP RHS
% 15 states = [12 aircraft; DELTAZ_act; DELTAW_act; PLA_act]
% ============================================================

function dx15 = rhs_trim_regulation_local(~, x15, x_trim12, u_trim, Kx, db, perf, ix_lqr)

x12 = x15(1:12);

xd10 = x12(ix_lqr) - x_trim12(ix_lqr);

du = -Kx * xd10;

u = u_trim;
u.DELTAZ = saturate_deg_local(u_trim.DELTAZ + du(1), perf.delta_max_deg);
u.DELTAW = saturate_deg_local(u_trim.DELTAW + du(2), perf.delta_max_deg);
u.PLA    = u_trim.PLA;

dx15 = ice_step_6dof(0, x15, u, db, perf);

end

% ============================================================
% RECOVER COMMAND SIGNALS
% ============================================================

function [delz_cmd, delw_cmd, pla_cmd] = recover_trim_regulation_signals( ...
    x15_hist, x_trim12, u_trim, Kx, perf, ix_lqr)

n = size(x15_hist,1);

delz_cmd = zeros(n,1);
delw_cmd = zeros(n,1);
pla_cmd  = zeros(n,1);

for i = 1:n
    x15 = x15_hist(i,1:15).';
    x12 = x15(1:12);

    xd10 = x12(ix_lqr) - x_trim12(ix_lqr);
    du = -Kx * xd10;

    delz_cmd(i) = saturate_deg_local(u_trim.DELTAZ + du(1), perf.delta_max_deg);
    delw_cmd(i) = saturate_deg_local(u_trim.DELTAW + du(2), perf.delta_max_deg);
    pla_cmd(i)  = u_trim.PLA;
end

end

% ============================================================
% GAMMA / AIRDATA
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

R_b2n = [ cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
          cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
             -sth,                sphi*cth,                 cphi*cth      ];

Vn = R_b2n * [u; v; w];
gamma_deg = rad2deg(atan2(-Vn(3), sqrt(Vn(1)^2 + Vn(2)^2)));

end

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

for kk = 2:numel(s)
    if s(kk) == 0
        s(kk) = s(kk-1);
    end
end

for kk = numel(s)-1:-1:1
    if s(kk) == 0
        s(kk) = s(kk+1);
    end
end

extrema_idx_local = [];
for kk = 1:numel(s)-1
    if s(kk) ~= 0 && s(kk+1) ~= 0 && s(kk+1) ~= s(kk)
        extrema_idx_local(end+1,1) = kk+1; %#ok<AGROW>
    end
end

if ~isempty(extrema_idx_local)
    keep = true(size(extrema_idx_local));
    for kk = 2:numel(extrema_idx_local)
        if extrema_idx_local(kk) - extrema_idx_local(kk-1) <= 1
            keep(kk) = false;
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

% ============================================================
% SIM HEALTH CHECK
% ============================================================

function [sim_ok, fail_reason] = sim_health_check_local(t, x15, tfinal, ode_fail, ode_fail_msg)

sim_ok = true;
fail_reason = 'ok';

if nargin >= 4 && ode_fail
    sim_ok = false;
    fail_reason = ['ode45_failed: ' ode_fail_msg];
    return
end

if isempty(t) || isempty(x15)
    sim_ok = false;
    fail_reason = 'empty_output';
    return
end

if any(~isfinite(t)) || any(~isfinite(x15(:)))
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
% SATURATION HELPER
% ============================================================

function y = saturate_deg_local(x, lim)
y = min(max(x, -lim), lim);
end

% ============================================================
% CSV EXPORT
% ============================================================

function save_case_csv_local(sim)

case_name_safe = lower(string(sim.name));
case_name_safe = regexprep(case_name_safe, '[^a-zA-Z0-9]+', '_');
case_name_safe = regexprep(case_name_safe, '^_+|_+$', '');

fname = sprintf('lqr_tvc_%s_vt_%g_alt_%g.csv', case_name_safe, sim.case.Vt, sim.case.alt_m);
fname = regexprep(fname, '[^a-zA-Z0-9._-]+', '_');

if ~isempty(sim.t_cl) && ~isempty(sim.x_cl)
    x = sim.x_cl;

    n = numel(sim.t_cl);
    trim_ok_col = repmat(sim.trim_ok, n, 1);
    sim_ok_col = repmat(sim.sim_ok, n, 1);
    conv_col = repmat(sim.is_converging, n, 1);

    T = table( ...
        sim.t_cl(:), ...
        trim_ok_col, ...
        sim_ok_col, ...
        conv_col, ...
        x(:,1), x(:,2), x(:,3), ...
        x(:,4), x(:,5), x(:,6), ...
        rad2deg(x(:,7)), rad2deg(x(:,8)), rad2deg(x(:,9)), ...
        x(:,10), x(:,11), x(:,12), ...
        sim.Vt_cl(:), sim.alpha_cl(:), sim.beta_cl(:), sim.gamma_cl(:), ...
        sim.alt_cl(:), ...
        sim.delz_cmd(:), sim.delw_cmd(:), sim.pla_cmd(:), ...
        sim.DELTAZ_act(:), sim.DELTAW_act(:), sim.PLA_act(:), ...
        'VariableNames', { ...
        't_s', ...
        'trim_ok','sim_ok','is_converging', ...
        'u_mps','v_mps','w_mps', ...
        'p_rad_s','q_rad_s','r_rad_s', ...
        'phi_deg','theta_deg','psi_deg', ...
        'xN_m','yE_m','zD_m', ...
        'Vt_mps','alpha_deg','beta_deg','gamma_deg', ...
        'altitude_m', ...
        'DELTAZ_cmd_deg','DELTAW_cmd_deg','PLA_cmd_deg', ...
        'DELTAZ_act_deg','DELTAW_act_deg','PLA_act_deg'});

    writetable(T, fname);
    fprintf('Saved %s\n', fname);
    return
end

T = table( ...
    string(sim.name), ...
    sim.case.Vt, ...
    sim.case.alt_m, ...
    sim.trim_ok, ...
    sim.sim_ok, ...
    string(sim.fail_reason), ...
    sim.trim_info.alpha_deg, ...
    sim.trim_info.theta_deg, ...
    sim.trim_info.gamma_deg, ...
    sim.trim_info.PLA_deg, ...
    sim.trim_info.DELTAZ_deg, ...
    sim.trim_info.resnorm, ...
    'VariableNames', { ...
    'case_name','Vt_mps','alt_m', ...
    'trim_ok','sim_ok','fail_reason', ...
    'alpha_deg','theta_deg','gamma_deg','PLA_deg','DELTAZ_deg','trim_resnorm'});

writetable(T, fname);
fprintf('Saved summary %s\n', fname);

end

% ============================================================
% PLOTTING
% ============================================================

function plot_case_results_local(sim)

if isempty(sim.t_cl) || isempty(sim.x_cl)
    return
end

figure('Color','w','Name',sim.name);
tiledlayout(4,2);

nexttile
plot(sim.t_cl, sim.Vt_cl, 'LineWidth',1.4);
grid on
ylabel('V_t (m/s)')
title(sim.name)

nexttile
plot(sim.t_cl, sim.alpha_cl, 'LineWidth',1.4);
grid on
ylabel('\alpha (deg)')

nexttile
plot(sim.t_cl, sim.alt_cl, 'LineWidth',1.4); hold on
yline(sim.alt_trim,'--');
if ~isempty(sim.extrema_idx)
    plot(sim.t_cl(sim.extrema_idx), sim.alt_cl(sim.extrema_idx), 'ro', 'MarkerSize', 6, 'LineWidth', 1.2);
end
grid on
ylabel('h (m)')

nexttile
plot(sim.t_cl, rad2deg(sim.x_cl(:,8)), 'LineWidth',1.4);
grid on
ylabel('\theta (deg)')

nexttile
plot(sim.t_cl, sim.delz_cmd, '--', 'LineWidth',1.2); hold on
plot(sim.t_cl, sim.DELTAZ_act, 'LineWidth',1.4);
grid on
ylabel('\Delta Z (deg)')
legend('\Delta Z_{cmd}','\Delta Z_{act}','Location','best')

nexttile
plot(sim.t_cl, sim.delw_cmd, '--', 'LineWidth',1.2); hold on
plot(sim.t_cl, sim.DELTAW_act, 'LineWidth',1.4);
grid on
ylabel('\Delta W (deg)')
legend('\Delta W_{cmd}','\Delta W_{act}','Location','best')

nexttile
plot(sim.t_cl, sim.PLA_act, 'LineWidth',1.4);
grid on
ylabel('PLA (deg)')
xlabel('t (s)')

nexttile
plot(sim.t_cl, sim.gamma_cl, 'LineWidth',1.4); hold on
plot(sim.t_cl, rad2deg(sim.x_cl(:,7)), 'LineWidth',1.2);
plot(sim.t_cl, rad2deg(sim.x_cl(:,9)), 'LineWidth',1.2);
grid on
ylabel('angle (deg)')
xlabel('t (s)')
legend('\gamma','\phi','\psi','Location','best')

end