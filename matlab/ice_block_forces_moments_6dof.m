function FM = ice_block_forces_moments_6dof(state, aero, eng, eff, db, perf)
% ICE_BLOCK_FORCES_MOMENTS_6DOF
% Combines aerodynamic and engine/TVC forces and moments in BODY axes.
%
% Inputs:
%   state   struct with rigid-body state
%   aero    struct with aerodynamic coefficients / qbar
%   eng     struct with engine outputs
%   eff     struct with control inputs:
%             eff.DELTAZ [deg]
%             eff.DELTAW [deg]
%             eff.PLA
%   db      database struct
%   perf    performance / actuator struct
%
% Outputs:
%   FM.Fb [3x1] body-axis total force [N]
%   FM.Mb [3x1] body-axis total moment [N*m]

if nargin < 6 || isempty(perf)
    perf = default_tvc_perf();
end

qbar = aero.qbar;
S = db.S;
b = db.b;
c = db.c;

% -------------------------------------------------
% Aerodynamic body-axis forces
% -------------------------------------------------
Faero_b = qbar * S * [aero.CXb; aero.CYb; aero.CZb];

% -------------------------------------------------
% Aerodynamic body-axis moments
% -------------------------------------------------
Maero_b = qbar * S * [b * aero.Cl;
                      c * aero.Cm;
                      b * aero.Cn];

% -------------------------------------------------
% Engine thrust in body axes with TVC effectiveness / loss
% -------------------------------------------------
T = eng.FG;

dz = deg2rad(eff.DELTAZ);
dw = deg2rad(eff.DELTAW);

delta_mag_deg = hypot(eff.DELTAZ, eff.DELTAW);

if isfield(perf, 'delta_loss_ref_deg')
    delta_ref_deg = perf.delta_loss_ref_deg;
else
    delta_ref_deg = perf.delta_max_deg;
end

eta_x = 1 - perf.k_loss * (delta_mag_deg / max(delta_ref_deg, 1e-6))^2;
eta_x = max(0.5, eta_x);

Tx =  T * eta_x        * cos(dz) * cos(dw);
Ty =  T * perf.eta_tvc * sin(dw);
Tz = -T * perf.eta_tvc * sin(dz) * cos(dw);

Fthrust_b = [Tx; Ty; Tz];

% -------------------------------------------------
% Engine drag in body axes
% -------------------------------------------------
Fdrag_eng_b = [-eng.DD; 0; 0];

% -------------------------------------------------
% Moments from engine force application points
% -------------------------------------------------
rFG = [db.xjFG; 0; db.zjFG];
rDD = [db.xjDD; 0; db.zjDD];

Mthrust_b = cross(rFG, Fthrust_b);
Mdrag_b   = cross(rDD, Fdrag_eng_b);

% -------------------------------------------------
% Total forces and moments
% -------------------------------------------------
Fb = Faero_b + Fthrust_b + Fdrag_eng_b;
Mb = Maero_b + Mthrust_b + Mdrag_b;

FM = struct();
FM.Fb = Fb;
FM.Mb = Mb;

FM.Faero_b   = Faero_b;
FM.Fthrust_b = Fthrust_b;
FM.Fdrag_b   = Fdrag_eng_b;

FM.Maero_b   = Maero_b;
FM.Mthrust_b = Mthrust_b;
FM.Mdrag_b   = Mdrag_b;

FM.eta_x   = eta_x;
FM.eta_tvc = perf.eta_tvc;

end