function sim = simulate_saved_turn_trim_lqr_6dof(csvFile, Vt_target, alt_target)
% SIMULATE_SAVED_TURN_TRIM_LQR_6DOF
%
% Reads one trimmed coordinated-turn point from a saved CSV, linearises at
% that point, designs an LQR using TVC only, and simulates the full
% nonlinear 15-state plant with actuator dynamics.
%
% IMPORTANT:
% - Regulates only states that should remain approximately fixed in a
%   steady coordinated turn:
%       [u v w p q r phi theta]
% - Does NOT regulate psi to a constant value, because psi should drift
%   during a steady turn.
% - PLA is held fixed at the trimmed value.
% - Saves a CSV of the time history automatically.
%
% Inputs:
%   csvFile    : CSV from tvc_turn_envelope_6dof... script
%   Vt_target  : target TAS [m/s] used to select nearest feasible row
%   alt_target : target altitude [m]
%
% Example:
%   sim = simulate_saved_turn_trim_lqr_6dof( ...
%       'tvc_turn_envelope_6dof_psidot_tas.csv', 220, 10000);

if nargin < 1 || isempty(csvFile)
    csvFile = 'tvc_turn_envelope_6dof_psidot_tas.csv';
end
if nargin < 2
    Vt_target = 402;
end
if nargin < 3
    alt_target = 10000;
end

db = ice_init();

% ============================================================
% TVC / engine dynamics performance model
% ============================================================
perf = default_tvc_perf();
perf.delta_max_deg = 25;
perf.delta_dot_max_deg_s = 50;
perf.tau_tvc_s = 0.1;
perf.tau_thrust_s = 0.42;
perf.PLA_dot_max_per_s = inf;

% ============================================================
% LQR on [u v w p q r phi theta]
% ============================================================
Q = diag([ ...
    1e-1 ... % u
    2e-1 ... % v
    2e-1 ... % w
    5e0  ... % p
    8e0  ... % q
    5e0  ... % r
    6e0  ... % phi
    6e0  ... % theta
]);

R = diag([ ...
    1e-1 ... % DELTAZ
    5e-1 ... % DELTAW
]);

ode_opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.02);

% ============================================================
% Read CSV and choose nearest feasible point
% ============================================================
T = readtable(csvFile);

requiredVars = { ...
    'Altitude_m','Feasible', ...
    'alpha_deg','beta_deg','theta_deg', ...
    'p_deg_s','q_deg_s','r_deg_s', ...
    'PLA_deg','DELTAZ_deg','DELTAW_deg'};

for k = 1:numel(requiredVars)
    if ~ismember(requiredVars{k}, T.Properties.VariableNames)
        error('Missing required column "%s" in %s', requiredVars{k}, csvFile);
    end
end

if ismember('Vt_mps', T.Properties.VariableNames)
    Vname = 'Vt_mps';
elseif ismember('TAS_mps', T.Properties.VariableNames)
    Vname = 'TAS_mps';
else
    error('CSV must contain either Vt_mps or TAS_mps.');
end

idx = T.Feasible == 1 & ...
      isfinite(T.(Vname)) & ...
      isfinite(T.Altitude_m);

if ~any(idx)
    error('No feasible rows found in %s', csvFile);
end

Vall = T.(Vname);
Hall = T.Altitude_m;

d2 = inf(height(T),1);
d2(idx) = (Vall(idx) - Vt_target).^2 + ((Hall(idx) - alt_target)/100).^2;
[~, iBest] = min(d2);

row = T(iBest,:);

fprintf('Selected row from %s\n', csvFile);
fprintf('  Target request : Vt = %.3f m/s, alt = %.3f m\n', Vt_target, alt_target);
fprintf('  Selected point : Vt = %.3f m/s, alt = %.3f m\n', row.(Vname), row.Altitude_m);

% ============================================================
% Build exact 12-state trim from saved row
% ============================================================
Vt = row.(Vname);
alt_m = row.Altitude_m;

alpha_rad = deg2rad(row.alpha_deg);
beta_rad  = deg2rad(row.beta_deg);
theta_rad = deg2rad(row.theta_deg);

phi_rad = 0;
if ismember('phi_deg', T.Properties.VariableNames) && isfinite(row.phi_deg)
    phi_rad = deg2rad(row.phi_deg);
end

p_rad_s = deg2rad(row.p_deg_s);
q_rad_s = deg2rad(row.q_deg_s);
r_rad_s = deg2rad(row.r_deg_s);

u_b = Vt * cos(alpha_rad) * cos(beta_rad);
v_b = Vt * sin(beta_rad);
w_b = Vt * sin(alpha_rad) * cos(beta_rad);

x_trim12 = [ ...
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
u_trim.PLA    = row.PLA_deg;
u_trim.DELTAZ = row.DELTAZ_deg;
u_trim.DELTAW = row.DELTAW_deg;

% ============================================================
% Build 15-state dynamic plant IC
% x15 = [12 rigid-body states; DELTAZ_act; DELTAW_act; PLA_act]
% ============================================================
x_trim15 = [x_trim12; u_trim.DELTAZ; u_trim.DELTAW; u_trim.PLA];

% ============================================================
% Check derivatives at trim
% ============================================================
dx12 = ice_step_6dof(0, x_trim12, u_trim, db);
dx15 = ice_step_6dof(0, x_trim15, u_trim, db, perf);

fprintf('12-state trim check norm = %.3e\n', norm(dx12(1:12)));
fprintf('15-state trim check norm = %.3e\n', norm(dx15(1:15)));

% ============================================================
% Linearise at turn trim
% ============================================================
[A,B] = linearize_ice_6dof(x_trim12, u_trim, db);

% Reduced regulation state: [u v w p q r phi theta]
xRegIdx = 1:8;

% Assumed input ordering in linearize_ice_6dof: [PLA DELTAZ DELTAW]
uTvcIdx = [2 3];

Ared = A(xRegIdx, xRegIdx);
Bred = B(xRegIdx, uTvcIdx);

K = lqr(Ared, Bred, Q, R);

fprintf('Closed-loop reduced eigenvalues for turn-hold inner loop:\n');
disp(eig(Ared - Bred*K));

% ============================================================
% Simulate with perturbation
% ============================================================
x0_12 = x_trim12;

% then you perturb it
x0_12(5) = x0_12(5) + deg2rad(5);   % q
x0_12(8) = x0_12(8) + deg2rad(5);   % theta
x0_12(7) = x0_12(7) + deg2rad(5);   % phi

x0_15 = [x0_12; u_trim.DELTAZ; u_trim.DELTAW; u_trim.PLA];

rhs = @(t,x) rhs_turn_trim_regulation_local(t, x, x_trim12, u_trim, K, db, perf);

fprintf('Starting nonlinear sim...\n');
[t, x15_hist] = ode45(rhs, [0 120], x0_15, ode_opts);

x12_hist = x15_hist(:,1:12);

DELTAZ_act = x15_hist(:,13);
DELTAW_act = x15_hist(:,14);
PLA_act    = x15_hist(:,15);

[Vt_hist, alpha_hist, beta_hist, gamma_hist] = derive_airdata_gamma_from_states_local(x12_hist);
phi_hist = rad2deg(x12_hist(:,7));
theta_hist = rad2deg(x12_hist(:,8));
psi_hist = rad2deg(x12_hist(:,9));
alt_hist = -x12_hist(:,12);

[delz_cmd, delw_cmd, pla_cmd] = recover_turn_regulation_signals_local( ...
    x15_hist, x_trim12, u_trim, K, perf);

% ============================================================
% Save outputs in struct
% ============================================================
sim = struct();
sim.csvFile = csvFile;
sim.selectedRow = row;
sim.x_trim12 = x_trim12;
sim.x_trim15 = x_trim15;
sim.u_trim = u_trim;
sim.dx12 = dx12;
sim.dx15 = dx15;
sim.A = A;
sim.B = B;
sim.K = K;
sim.t = t;
sim.x15_hist = x15_hist;
sim.x12_hist = x12_hist;
sim.Vt_hist = Vt_hist;
sim.alpha_hist = alpha_hist;
sim.beta_hist = beta_hist;
sim.gamma_hist = gamma_hist;
sim.phi_hist = phi_hist;
sim.theta_hist = theta_hist;
sim.psi_hist = psi_hist;
sim.alt_hist = alt_hist;
sim.DELTAZ_act = DELTAZ_act;
sim.DELTAW_act = DELTAW_act;
sim.PLA_act = PLA_act;
sim.DELTAZ_cmd = delz_cmd;
sim.DELTAW_cmd = delw_cmd;
sim.PLA_cmd = pla_cmd;

% ============================================================
% Plots
% ============================================================
figure('Color','w');
tiledlayout(4,2)

nexttile
plot(t, Vt_hist, 'LineWidth', 1.4)
grid on
ylabel('V_t (m/s)')
title(sprintf('LQR turn trim hold | V_t = %.1f m/s | h = %.0f m', Vt, alt_m))

nexttile
plot(t, gamma_hist, 'LineWidth', 1.4)
grid on
ylabel('\gamma (deg)')

nexttile
plot(t, phi_hist, 'LineWidth', 1.4)
grid on
ylabel('\phi (deg)')

nexttile
plot(t, psi_hist, 'LineWidth', 1.4)
grid on
ylabel('\psi (deg)')

nexttile
plot(t, theta_hist, 'LineWidth', 1.4)
grid on
ylabel('\theta (deg)')

nexttile
plot(t, alpha_hist, 'LineWidth', 1.4)
grid on
ylabel('\alpha (deg)')

nexttile
plot(t, delz_cmd - u_trim.DELTAZ, '--', 'LineWidth', 1.2); hold on
plot(t, DELTAZ_act - u_trim.DELTAZ, 'LineWidth', 1.4)
plot(t, delw_cmd - u_trim.DELTAW, ':', 'LineWidth', 1.2)
plot(t, DELTAW_act - u_trim.DELTAW, '-.', 'LineWidth', 1.4)
grid on
ylabel('TVC dev. (deg)')
legend('\Delta Z_{cmd}-\Delta Z_{trim}', ...
       '\Delta Z_{act}-\Delta Z_{trim}', ...
       '\Delta W_{cmd}-\Delta W_{trim}', ...
       '\Delta W_{act}-\Delta W_{trim}', ...
       'Location','best')
xlabel('Time (s)')

nexttile
plot(t, alt_hist, 'LineWidth', 1.4)
grid on
ylabel('Altitude (m)')
xlabel('Time (s)')

% ============================================================
% SAVE CSV
% ============================================================
fname = sprintf('lqr_turn_vt_%.0f_alt_%.0f.csv', Vt, alt_m);

Tsave = table( ...
    t(:), ...
    Vt_hist(:), ...
    alpha_hist(:), ...
    beta_hist(:), ...
    gamma_hist(:), ...
    phi_hist(:), ...
    theta_hist(:), ...
    psi_hist(:), ...
    alt_hist(:), ...
    DELTAZ_act(:), ...
    DELTAW_act(:), ...
    PLA_act(:), ...
    delz_cmd(:), ...
    delw_cmd(:), ...
    pla_cmd(:), ...
    DELTAZ_act(:) - u_trim.DELTAZ, ...
    DELTAW_act(:) - u_trim.DELTAW, ...
    delz_cmd(:) - u_trim.DELTAZ, ...
    delw_cmd(:) - u_trim.DELTAW, ...
    'VariableNames', { ...
    't_s', ...
    'Vt_mps', ...
    'alpha_deg', ...
    'beta_deg', ...
    'gamma_deg', ...
    'phi_deg', ...
    'theta_deg', ...
    'psi_deg', ...
    'altitude_m', ...
    'DELTAZ_act_deg', ...
    'DELTAW_act_deg', ...
    'PLA_act_deg', ...
    'DELTAZ_cmd_deg', ...
    'DELTAW_cmd_deg', ...
    'PLA_cmd_deg', ...
    'DELTAZ_act_dev_deg', ...
    'DELTAW_act_dev_deg', ...
    'DELTAZ_cmd_dev_deg', ...
    'DELTAW_cmd_dev_deg'});

writetable(Tsave, fname);
fprintf('Saved CSV: %s\n', fname);

end

% =====================================================================
% CLOSED-LOOP RHS FOR TURN TRIM HOLD
% 15 states = [12 rigid-body; DELTAZ_act; DELTAW_act; PLA_act]
% =====================================================================

function dx15 = rhs_turn_trim_regulation_local(~, x15, x_trim12, u_trim, K, db, perf)

x12 = x15(1:12);

% Regulate only states that should remain approximately fixed in a steady turn
% [u v w p q r phi theta]
xd = x12(1:8) - x_trim12(1:8);

du = -K * xd;

u = u_trim;
u.DELTAZ = saturate_deg_local(u_trim.DELTAZ + du(1), perf.delta_max_deg);
u.DELTAW = saturate_deg_local(u_trim.DELTAW + du(2), perf.delta_max_deg);

% Keep PLA fixed
u.PLA = u_trim.PLA;

dx15 = ice_step_6dof(0, x15, u, db, perf);

end

% =====================================================================
% RECOVER COMMAND SIGNALS
% =====================================================================

function [delz_cmd, delw_cmd, pla_cmd] = recover_turn_regulation_signals_local( ...
    x15_hist, x_trim12, u_trim, K, perf)

n = size(x15_hist,1);

delz_cmd = zeros(n,1);
delw_cmd = zeros(n,1);
pla_cmd = zeros(n,1);

for i = 1:n
    x12 = x15_hist(i,1:12).';
    xd = x12(1:8) - x_trim12(1:8);
    du = -K * xd;

    delz_cmd(i) = saturate_deg_local(u_trim.DELTAZ + du(1), perf.delta_max_deg);
    delw_cmd(i) = saturate_deg_local(u_trim.DELTAW + du(2), perf.delta_max_deg);
    pla_cmd(i) = u_trim.PLA;
end

end

% =====================================================================
% HELPERS
% =====================================================================

function y = saturate_deg_local(x, lim)
y = min(max(x, -lim), lim);
end

function [Vt, alpha_deg, beta_deg, gamma_deg] = derive_airdata_gamma_from_states_local(xhist)

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