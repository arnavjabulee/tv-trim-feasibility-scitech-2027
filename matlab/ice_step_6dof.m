function xdot = ice_step_6dof(~, x, u, db, perf)
% ICE_STEP_6DOF
%
% State order (new):
% [u v w p q r phi theta psi xN yE zD DELTAZ_act DELTAW_act PLA_act]
%
% Backward compatible:
% If x has only 12 states, old static behaviour is used:
% [u v w p q r phi theta psi xN yE zD]

if nargin < 5 || isempty(perf)
    perf = default_tvc_perf();
end

% ----------------------------
% Rigid-body states
% ----------------------------
state = struct();
state.u = x(1);
state.v = x(2);
state.w = x(3);

state.p = x(4);
state.q = x(5);
state.r = x(6);

state.phi   = x(7);
state.theta = x(8);
state.psi   = x(9);

state.xN = x(10);
state.yE = x(11);
state.zD = x(12);

Vt = norm(x(1:3));

% ----------------------------
% Backward-compatible static mode
% ----------------------------
if numel(x) < 15
    eng = ice_block_engine(state.zD, Vt, u.PLA, db);
    state.rho = eng.rho;

    aero = ice_block_aero_data_6dof(state, db);
    FM   = ice_block_forces_moments_6dof(state, aero, eng, u, db, perf);

    xdot = ice_block_eom_6dof(state, FM, db);
    return
end

% ----------------------------
% Dynamic actuator / throttle states
% ----------------------------
deltaZ_act = x(13);
deltaW_act = x(14);
PLA_act    = x(15);

% ----------------------------
% Command saturation
% ----------------------------
deltaZ_cmd = min(max(u.DELTAZ, -perf.delta_max_deg), perf.delta_max_deg);
deltaW_cmd = min(max(u.DELTAW, -perf.delta_max_deg), perf.delta_max_deg);

PLA_min = min(db.eng.PLA);
PLA_max = max(db.eng.PLA);
PLA_cmd = min(max(u.PLA, PLA_min), PLA_max);

% ----------------------------
% First-order desired derivatives
% ----------------------------
tau_tvc = max(perf.tau_tvc_s, 1e-6);
tau_thr = max(perf.tau_thrust_s, 1e-6);

deltaZ_dot_des = (deltaZ_cmd - deltaZ_act) / tau_tvc;
deltaW_dot_des = (deltaW_cmd - deltaW_act) / tau_tvc;
PLA_dot_des    = (PLA_cmd    - PLA_act   ) / tau_thr;

% ----------------------------
% Rate-limited derivatives
% Continuous-time limit for ODE use
% ----------------------------
deltaZ_dot = min(max(deltaZ_dot_des, -perf.delta_dot_max_deg_s), perf.delta_dot_max_deg_s);
deltaW_dot = min(max(deltaW_dot_des, -perf.delta_dot_max_deg_s), perf.delta_dot_max_deg_s);

if isfinite(perf.PLA_dot_max_per_s)
    PLA_dot = min(max(PLA_dot_des, -perf.PLA_dot_max_per_s), perf.PLA_dot_max_per_s);
else
    PLA_dot = PLA_dot_des;
end

% ----------------------------
% Use ACTUAL lagged states for engine and TVC
% ----------------------------
eng = ice_block_engine(state.zD, Vt, PLA_act, db);
state.rho = eng.rho;

aero = ice_block_aero_data_6dof(state, db);

eff_dyn = u;
eff_dyn.DELTAZ = deltaZ_act;
eff_dyn.DELTAW = deltaW_act;
eff_dyn.PLA    = PLA_act;

FM = ice_block_forces_moments_6dof(state, aero, eng, eff_dyn, db, perf);

% ----------------------------
% Rigid-body dynamics + actuator dynamics
% ----------------------------
xdot_rb = ice_block_eom_6dof(state, FM, db);

xdot = [xdot_rb;
        deltaZ_dot;
        deltaW_dot;
        PLA_dot];
end