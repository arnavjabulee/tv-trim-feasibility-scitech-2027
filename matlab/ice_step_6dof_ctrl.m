function dx = ice_step_6dof_ctrl(~, x, ucmd, db)
% ICE_STEP_6DOF_CTRL
%
% 6-DOF equations of motion using aerodynamic controls only.
%
% State:
%   x = [u v w p q r phi theta psi xN yE zD]'
%
% Control struct fields used:
%   PLA, PFLAP, SPLRL,
%   IBLFL, OBLFL, ELEVL, CLAML, AMTL, DEFL,
%   IBLFR, OBLFR, ELEVR, SPLRR, CLAMR, AMTR, DEFR
%
% TVC is forced to zero in this branch.

% -----------------------------
% Unpack state
% -----------------------------
u = x(1); v = x(2); w = x(3);
p = x(4); q = x(5); r = x(6);
phi = x(7); theta = x(8); psi = x(9);
zD = x(12);

alt_m = -zD;

% -----------------------------
% Atmosphere
% -----------------------------
rho = interp1(db.atm.h_m, db.atm.rho, alt_m, 'linear', 'extrap');
T_atm = interp1(db.atm.h_m, db.atm.T, alt_m, 'linear', 'extrap');

state = struct();
state.u   = u;
state.v   = v;
state.w   = w;
state.p   = p;
state.q   = q;
state.r   = r;
state.rho = rho;

% -----------------------------
% Controls with safe defaults
% -----------------------------
eff = struct();

eff.PLA   = get_ctrl(ucmd, 'PLA', 0);

eff.PFLAP = get_ctrl(ucmd, 'PFLAP', 0);
eff.SPLRL = get_ctrl(ucmd, 'SPLRL', 0);

eff.IBLFL = get_ctrl(ucmd, 'IBLFL', 0);
eff.OBLFL = get_ctrl(ucmd, 'OBLFL', 0);
eff.ELEVL = get_ctrl(ucmd, 'ELEVL', 0);
eff.CLAML = get_ctrl(ucmd, 'CLAML', 0);
eff.AMTL  = get_ctrl(ucmd, 'AMTL',  0);
eff.DEFL  = get_ctrl(ucmd, 'DEFL',  0);

eff.IBLFR = get_ctrl(ucmd, 'IBLFR', 0);
eff.OBLFR = get_ctrl(ucmd, 'OBLFR', 0);
eff.ELEVR = get_ctrl(ucmd, 'ELEVR', 0);
eff.SPLRR = get_ctrl(ucmd, 'SPLRR', 0);
eff.CLAMR = get_ctrl(ucmd, 'CLAMR', 0);
eff.AMTR  = get_ctrl(ucmd, 'AMTR',  0);
eff.DEFR  = get_ctrl(ucmd, 'DEFR',  0);

% Force TVC off in this branch
eff.DELTAZ = 0;
eff.DELTAW = 0;

% -----------------------------
% Aerodynamics
% -----------------------------
aero = ice_block_aero_data_6dof_ctrl(state, eff, db);

qbar = aero.qbar;

Fax = qbar * db.S * aero.CXb;
Fay = qbar * db.S * aero.CYb;
Faz = qbar * db.S * aero.CZb;

L = qbar * db.S * db.b * aero.Cl;
M = qbar * db.S * db.c * aero.Cm;
N = qbar * db.S * db.b * aero.Cn;

% -----------------------------
% Engine thrust / inlet drag
% -----------------------------
Vt = max(aero.Vt, 1e-6);
eng = ice_block_engine(zD, Vt, eff.PLA, db);

FG = eng.FG;
DD = eng.DD;

% Aero-controls branch: no thrust vectoring
T_body = [FG; 0; 0];
D_inlet_body = [-DD; 0; 0];

% -----------------------------
% Gravity in body axes
% -----------------------------
g = db.g;
m = db.m;

Xg = -m * g * sin(theta);
Yg =  m * g * sin(phi) * cos(theta);
Zg =  m * g * cos(phi) * cos(theta);

% -----------------------------
% Total forces in body axes
% -----------------------------
X = Fax + T_body(1) + D_inlet_body(1) + Xg;
Y = Fay + T_body(2) + D_inlet_body(2) + Yg;
Z = Faz + T_body(3) + D_inlet_body(3) + Zg;

% -----------------------------
% Translational equations
% -----------------------------
udot = r*v - q*w + X/m;
vdot = p*w - r*u + Y/m;
wdot = q*u - p*v + Z/m;

% -----------------------------
% Rotational equations
% -----------------------------
Ixx = db.Ixx; Iyy = db.Iyy; Izz = db.Izz; Ixz = db.Ixz;

Gamma  = Ixx*Izz - Ixz^2;
Gamma1 = Ixz*(Ixx - Iyy + Izz) / Gamma;
Gamma2 = (Izz*(Izz - Iyy) + Ixz^2) / Gamma;
Gamma3 = Izz / Gamma;
Gamma4 = Ixz / Gamma;
Gamma5 = (Izz - Ixx) / Iyy;
Gamma6 = Ixz / Iyy;
Gamma7 = ((Ixx - Iyy)*Ixx + Ixz^2) / Gamma;
Gamma8 = Ixx / Gamma;

pdot = Gamma1*p*q - Gamma2*q*r + Gamma3*L + Gamma4*N;
qdot = Gamma5*p*r - Gamma6*(p^2 - r^2) + M/Iyy;
rdot = Gamma7*p*q - Gamma1*q*r + Gamma4*L + Gamma8*N;

% -----------------------------
% Kinematics
% -----------------------------
cphi = cos(phi);   sphi = sin(phi);
cth  = cos(theta); sth  = sin(theta);
cpsi = cos(psi);   spsi = sin(psi);

% Guard against singularity
cth_safe = cth;
if abs(cth_safe) < 1e-6
    cth_safe = sign(cth_safe + eps) * 1e-6;
end
tth = sth / cth_safe;

phidot   = p + q*sphi*tth + r*cphi*tth;
thetadot = q*cphi - r*sphi;
psidot   = q*sphi/cth_safe + r*cphi/cth_safe;

R_b2n = [ cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
          cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
             -sth,                sphi*cth,                 cphi*cth      ];

Vn = R_b2n * [u; v; w];
xNdot = Vn(1);
yEdot = Vn(2);
zDdot = Vn(3);

dx = [ ...
    udot;
    vdot;
    wdot;
    pdot;
    qdot;
    rdot;
    phidot;
    thetadot;
    psidot;
    xNdot;
    yEdot;
    zDdot];
end

function val = get_ctrl(s, name, defaultVal)
if isfield(s, name) && ~isempty(s.(name)) && isfinite(s.(name))
    val = s.(name);
else
    val = defaultVal;
end
end