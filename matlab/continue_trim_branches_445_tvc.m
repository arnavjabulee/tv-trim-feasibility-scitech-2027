function branches = continue_trim_branches_445_tvc()
% CONTINUE_TRIM_BRANCHES_445_TVC
%
% Tracks the two known trim branches at Vt = 445 m/s for the TVC-only case:
%   - main-envelope branch
%   - upper-island branch
%
% This uses continuation in altitude:
%   solve one seed point in each branch
%   then march upward and downward using the previous converged solution
%   as the next initial guess
%
% This is meant to answer:
%   Are the "main envelope" and "island" actually separate trim branches?
%
% Assumptions:
%   - unrestricted PLA for feasibility
%   - DELTAZ still limited
%   - level trim only
%
% Unknowns:
%   y = [alpha_deg; theta_rad; PLA_deg; DELTAZ_deg]

db = ice_init();

% ============================================================
% USER SETTINGS
% ============================================================
Vt = 445;                     % m/s
alt_list_m = 0:100:20000;     % continuation grid
DELTAZ_lim = 35;              % deg
trim_tol = 1e-5;

showConsoleProgress = true;

% These two altitudes should lie inside the two feasible bands you already found.
seed_main_alt_m   = 12500;    % pick inside main-envelope band
seed_island_alt_m = 17500;    % pick inside island band

% Optional manual seed vectors:
% y = [alpha_deg; theta_rad; PLA_deg; DELTAZ_deg]
% Leave empty to use default initial guess at that altitude.
seed_main_y0   = [];
seed_island_y0 = [];

% Output CSV filenames
main_csv_name   = 'branch_445_main.csv';
island_csv_name = 'branch_445_island.csv';

% Initial-guess helper limits only
PLA_guess_min = min(db.eng.PLA);
PLA_guess_max = max(db.eng.PLA);

gamma_target_rad = 0.0;       % level trim

opts = optimoptions('fsolve', ...
    'Display','off', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-8, ...
    'MaxFunctionEvaluations',1200, ...
    'MaxIterations',150);

% ============================================================
% FIND SEED SOLUTIONS
% ============================================================
fprintf('\n============================================================\n');
fprintf('BRANCH CONTINUATION AT Vt = %.1f m/s\n', Vt);
fprintf('DELTAZ limit = %.1f deg | PLA unrestricted\n', DELTAZ_lim);
fprintf('============================================================\n');

if isempty(seed_main_y0)
    y0_main = default_trim_guess_level(Vt, seed_main_alt_m, PLA_guess_min, PLA_guess_max);
else
    y0_main = seed_main_y0(:);
end

if isempty(seed_island_y0)
    y0_island = default_trim_guess_level(Vt, seed_island_alt_m, PLA_guess_min, PLA_guess_max);
else
    y0_island = seed_island_y0(:);
end

[seed_main_sol, seed_main_info] = solve_trim_point_6dof( ...
    Vt, seed_main_alt_m, gamma_target_rad, db, y0_main, opts, trim_tol);

[seed_main_feas, ~, seed_main_reason] = classify_trim_result_tvc_unrestricted_pla( ...
    seed_main_info, seed_main_sol(4), DELTAZ_lim, trim_tol);

if ~seed_main_feas
    error('Main-branch seed at alt = %.1f m failed: %s', seed_main_alt_m, seed_main_reason);
end

[seed_island_sol, seed_island_info] = solve_trim_point_6dof( ...
    Vt, seed_island_alt_m, gamma_target_rad, db, y0_island, opts, trim_tol);

[seed_island_feas, ~, seed_island_reason] = classify_trim_result_tvc_unrestricted_pla( ...
    seed_island_info, seed_island_sol(4), DELTAZ_lim, trim_tol);

if ~seed_island_feas
    error('Island-branch seed at alt = %.1f m failed: %s', seed_island_alt_m, seed_island_reason);
end

fprintf('\nSeed solutions found.\n');
fprintf('  Main seed   at %.1f m: alpha=%.3f, theta=%.3f deg, PLA=%.3f, DELTAZ=%.3f\n', ...
    seed_main_alt_m, seed_main_sol(1), rad2deg(seed_main_sol(2)), seed_main_sol(3), seed_main_sol(4));
fprintf('  Island seed at %.1f m: alpha=%.3f, theta=%.3f deg, PLA=%.3f, DELTAZ=%.3f\n', ...
    seed_island_alt_m, seed_island_sol(1), rad2deg(seed_island_sol(2)), seed_island_sol(3), seed_island_sol(4));

% ============================================================
% CONTINUE BOTH BRANCHES
% ============================================================
main_branch = continue_one_branch( ...
    db, Vt, alt_list_m, gamma_target_rad, DELTAZ_lim, trim_tol, opts, ...
    seed_main_alt_m, seed_main_sol, 'main', showConsoleProgress);

island_branch = continue_one_branch( ...
    db, Vt, alt_list_m, gamma_target_rad, DELTAZ_lim, trim_tol, opts, ...
    seed_island_alt_m, seed_island_sol, 'island', showConsoleProgress);

branches = struct();
branches.Vt = Vt;
branches.DELTAZ_lim_deg = DELTAZ_lim;
branches.main = main_branch;
branches.island = island_branch;

% ============================================================
% SAVE CSV FILES
% ============================================================
save_branch_csv(main_branch, main_csv_name);
save_branch_csv(island_branch, island_csv_name);

fprintf('\nSaved %s\n', main_csv_name);
fprintf('Saved %s\n', island_csv_name);

% ============================================================
% PLOTS
% ============================================================
plot_branches(branches);

% ============================================================
% COMMAND WINDOW SUMMARY
% ============================================================
fprintf('\n============================================================\n');
fprintf('BRANCH SUMMARY\n');
fprintf('============================================================\n');
print_branch_summary(main_branch,   'MAIN BRANCH');
print_branch_summary(island_branch, 'ISLAND BRANCH');

end

% =====================================================================
% CONTINUE ONE BRANCH FROM A SEED
% =====================================================================

function branch = continue_one_branch( ...
    db, Vt, alt_list_m, gamma_target_rad, DELTAZ_lim, trim_tol, opts, ...
    seed_alt_m, seed_sol, branch_name, showConsoleProgress)

nA = numel(alt_list_m);

alpha_req  = NaN(nA,1);
theta_req  = NaN(nA,1);
gamma_req  = NaN(nA,1);
PLA_req    = NaN(nA,1);
DELTAZ_req = NaN(nA,1);
T_req      = NaN(nA,1);

feasible   = false(nA,1);
evaluated  = false(nA,1);
failure_code   = zeros(nA,1);
failure_reason = strings(nA,1);
failure_reason(:) = "Not evaluated";

trim_info  = cell(nA,1);
x_trim_cell = cell(nA,1);
u_trim_cell = cell(nA,1);

[~, seed_idx] = min(abs(alt_list_m - seed_alt_m));

% ---------- store seed ----------
store_solution(seed_idx, seed_sol);

% ---------- upward continuation ----------
prev_sol = seed_sol;
for ia = seed_idx+1:nA
    alt_m = alt_list_m(ia);
    [ok, sol, info, failCode, failReason] = try_one_point( ...
        Vt, alt_m, gamma_target_rad, db, prev_sol, opts, trim_tol, DELTAZ_lim);

    if ok
        store_solution(ia, sol, info);
        prev_sol = sol;
    else
        store_failure(ia, failCode, failReason, info);
        if showConsoleProgress
            fprintf('[%s][UP] alt=%6.0f | fail | %s | stopping upward continuation\n', ...
                upper(branch_name), alt_m, failReason);
        end
        break
    end
end

% ---------- downward continuation ----------
prev_sol = seed_sol;
for ia = seed_idx-1:-1:1
    alt_m = alt_list_m(ia);
    [ok, sol, info, failCode, failReason] = try_one_point( ...
        Vt, alt_m, gamma_target_rad, db, prev_sol, opts, trim_tol, DELTAZ_lim);

    if ok
        store_solution(ia, sol, info);
        prev_sol = sol;
    else
        store_failure(ia, failCode, failReason, info);
        if showConsoleProgress
            fprintf('[%s][DN] alt=%6.0f | fail | %s | stopping downward continuation\n', ...
                upper(branch_name), alt_m, failReason);
        end
        break
    end
end

branch = struct();
branch.name = string(branch_name);
branch.seed_alt_m = seed_alt_m;
branch.alt_list_m = alt_list_m(:);
branch.alpha_req = alpha_req;
branch.theta_req = theta_req;
branch.gamma_req = gamma_req;
branch.PLA_req = PLA_req;
branch.DELTAZ_req = DELTAZ_req;
branch.T_req = T_req;
branch.feasible = feasible;
branch.evaluated = evaluated;
branch.failure_code = failure_code;
branch.failure_reason = failure_reason;
branch.trim_info = trim_info;
branch.x_trim_cell = x_trim_cell;
branch.u_trim_cell = u_trim_cell;

    function store_solution(ia, sol, info)
        if nargin < 3
            info = [];
        end

        if isempty(info)
            [~, info] = solve_trim_point_6dof( ...
                Vt, alt_list_m(ia), gamma_target_rad, db, sol, opts, trim_tol);
        end

        alpha_req(ia)  = sol(1);
        theta_req(ia)  = rad2deg(sol(2));
        gamma_req(ia)  = info.gamma_deg;
        PLA_req(ia)    = sol(3);
        DELTAZ_req(ia) = sol(4);
        T_req(ia)      = estimate_thrust_from_pla_N(db, alt_list_m(ia), Vt, sol(3));

        feasible(ia) = true;
        evaluated(ia) = true;
        failure_code(ia) = 1;
        failure_reason(ia) = "Feasible";

        trim_info{ia} = info;
        x_trim_cell{ia} = info.x_trim;
        u_trim_cell{ia} = info.u_trim;

        if showConsoleProgress
            fprintf('[%s] alt=%6.0f | feas=1 | alpha=%7.3f | theta=%7.3f | PLA=%9.3f | DELTAZ=%7.3f\n', ...
                upper(branch_name), alt_list_m(ia), sol(1), rad2deg(sol(2)), sol(3), sol(4));
        end
    end

    function store_failure(ia, failCode, failReason, info)
        feasible(ia) = false;
        evaluated(ia) = true;
        failure_code(ia) = failCode;
        failure_reason(ia) = failReason;
        trim_info{ia} = info;
    end
end

% =====================================================================
% SAVE BRANCH TO CSV
% =====================================================================

function save_branch_csv(branch, fname)

T = table( ...
    branch.alt_list_m(:), ...
    branch.alpha_req(:), ...
    branch.theta_req(:), ...
    branch.gamma_req(:), ...
    branch.PLA_req(:), ...
    branch.DELTAZ_req(:), ...
    branch.T_req(:), ...
    double(branch.feasible(:)), ...
    double(branch.evaluated(:)), ...
    branch.failure_code(:), ...
    cellstr(branch.failure_reason(:)), ...
    'VariableNames', { ...
    'Altitude_m', ...
    'Alpha_deg', ...
    'Theta_deg', ...
    'Gamma_deg', ...
    'PLA_deg', ...
    'DELTAZ_deg', ...
    'Treq_N', ...
    'Feasible', ...
    'Evaluated', ...
    'FailureCode', ...
    'FailureReason'});

writetable(T, fname);

end

% =====================================================================
% TRY ONE POINT
% =====================================================================

function [ok, sol, info, failCode, failReason] = try_one_point( ...
    Vt, alt_m, gamma_target_rad, db, y0, opts, trim_tol, DELTAZ_lim)

ok = false;
sol = nan(4,1);
info = struct();

try
    [sol, info] = solve_trim_point_6dof( ...
        Vt, alt_m, gamma_target_rad, db, y0, opts, trim_tol);

    [isFeasible, failCode, failReason] = classify_trim_result_tvc_unrestricted_pla( ...
        info, sol(4), DELTAZ_lim, trim_tol);

    ok = isFeasible;

catch ME
    failCode = 2;
    failReason = "Solver/lookup error";
    info = struct('error_message', ME.message);
end
end

% =====================================================================
% PLOTTING
% =====================================================================

function plot_branches(branches)

m = branches.main;
u = branches.island;

alt_m = m.alt_list_m;

figure('Color','w', 'Name', sprintf('Branch continuation at Vt = %.1f m/s', branches.Vt));
tiledlayout(2,3, 'TileSpacing','compact', 'Padding','compact');

nexttile;
hold on;
plot(m.PLA_req, alt_m, 'b-', 'LineWidth', 1.5);
plot(u.PLA_req, alt_m, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('PLA (deg)');
ylabel('Altitude (m)');
title('PLA');
legend('Main branch','Island branch','Location','best');

nexttile;
hold on;
plot(m.DELTAZ_req, alt_m, 'b-', 'LineWidth', 1.5);
plot(u.DELTAZ_req, alt_m, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('\Delta Z (deg)');
ylabel('Altitude (m)');
title('\Delta Z');

nexttile;
hold on;
plot(m.alpha_req, alt_m, 'b-',  'LineWidth', 1.5);
plot(m.theta_req, alt_m, 'b--', 'LineWidth', 1.2);
plot(u.alpha_req, alt_m, 'r-',  'LineWidth', 1.5);
plot(u.theta_req, alt_m, 'r--', 'LineWidth', 1.2);
grid on;
xlabel('Angle (deg)');
ylabel('Altitude (m)');
title('\alpha and \theta');
legend('\alpha main','\theta main','\alpha island','\theta island','Location','best');

nexttile;
hold on;
scatter(zeros(nnz(m.feasible),1) + 1.00, alt_m(m.feasible), 18, 'b', 'filled');
scatter(zeros(nnz(u.feasible),1) + 1.03, alt_m(u.feasible), 18, 'r', 'filled');
xlim([0.9 1.12]);
grid on;
xlabel('Feasible points');
title('Feasible branch segments');
ylabel('Altitude (m)');
legend('Main branch','Island branch','Location','best');

nexttile;
hold on;
plot(m.failure_code, alt_m, 'b.-', 'LineWidth', 1.2, 'MarkerSize', 12);
plot(u.failure_code, alt_m, 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12);
grid on;
xlabel('Failure code');
ylabel('Altitude (m)');
title('Failure code');
legend('Main branch','Island branch','Location','best');

nexttile;
hold on;
if all(isnan(m.T_req)) && all(isnan(u.T_req))
    text(0.12, 0.5, 'T_{req} not implemented', 'Units','normalized');
    axis off;
else
    plot(m.T_req, alt_m, 'b-', 'LineWidth', 1.5);
    plot(u.T_req, alt_m, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('T_{req} (N)');
    ylabel('Altitude (m)');
    title('Estimated thrust');
    legend('Main branch','Island branch','Location','best');
end

sgtitle(sprintf('Branch continuation at V_t = %.1f m/s, DELTAZ limit = %.1f deg', ...
    branches.Vt, branches.DELTAZ_lim_deg));

end

% =====================================================================
% SUMMARY
% =====================================================================

function print_branch_summary(branch, label)

idx = branch.feasible;

fprintf('\n%s\n', label);
fprintf('  feasible altitudes found: %d\n', nnz(idx));

if ~any(idx)
    fprintf('  no feasible points traced\n');
    return
end

alts = branch.alt_list_m(idx);
fprintf('  min altitude: %.1f m\n', min(alts));
fprintf('  max altitude: %.1f m\n', max(alts));

fprintf('  PLA range   : %.3f to %.3f deg\n', ...
    min(branch.PLA_req(idx)), max(branch.PLA_req(idx)));
fprintf('  DELTAZ range: %.3f to %.3f deg\n', ...
    min(branch.DELTAZ_req(idx)), max(branch.DELTAZ_req(idx)));
fprintf('  alpha range : %.3f to %.3f deg\n', ...
    min(branch.alpha_req(idx)), max(branch.alpha_req(idx)));
fprintf('  theta range : %.3f to %.3f deg\n', ...
    min(branch.theta_req(idx)), max(branch.theta_req(idx)));

end

% =====================================================================
% CORE SOLVER
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

alpha_deg  = y(1);
theta_rad  = y(2);
PLA_deg    = y(3);
DELTAZ_deg = y(4);

if abs(alpha_deg) > 40 || ...
   abs(rad2deg(theta_rad)) > 50 || ...
   abs(PLA_deg) > 1e4 || ...
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
% INITIAL GUESS
% =====================================================================

function y0 = default_trim_guess_level(Vt, alt_m, PLA_guess_min, PLA_guess_max)

alpha0_deg = max(1, min(10, 2 + 0.015*(200 - Vt) + 0.0001*alt_m));
theta0_rad = deg2rad(alpha0_deg);
PLA0_deg   = max(PLA_guess_min + 5, min(PLA_guess_max - 5, 50 + 0.004*alt_m + 0.04*(Vt - 100)));
DELTAZ0_deg = 0;

y0 = [alpha0_deg; theta0_rad; PLA0_deg; DELTAZ0_deg];

end

% =====================================================================
% GAMMA
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
% FEASIBILITY CLASSIFICATION
% =====================================================================

function [isFeasible, code, reason] = classify_trim_result_tvc_unrestricted_pla( ...
    info, DELTAZ_deg, DELTAZ_lim, trim_tol)

% 1 feasible
% 2 solver/lookup error
% 3 exitflag/non-convergence
% 4 DELTAZ limit
% 5 residual norm tolerance
% 6 other

isFeasible = false;

if ~(isfinite(DELTAZ_deg) && abs(DELTAZ_deg) <= DELTAZ_lim)
    code = 4;
    reason = "DELTAZ limit";
    return
end

if info.exitflag <= 0
    code = 3;
    reason = "Non-converged exitflag";
    return
end

if ~(isfinite(info.resnorm) && info.resnorm < trim_tol)
    code = 5;
    reason = "Residual norm tolerance";
    return
end

if info.success
    isFeasible = true;
    code = 1;
    reason = "Feasible";
    return
end

code = 6;
reason = "Other failure";
end

% =====================================================================
% OPTIONAL THRUST EXTRACTION HOOK
% =====================================================================

function Treq_N = estimate_thrust_from_pla_N(~, ~, ~, ~)
Treq_N = NaN;
end