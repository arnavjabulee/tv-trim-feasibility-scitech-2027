function aero = ice_block_aero_data_6dof(state, db)
% ICE_BLOCK_AERO_DATA_6DOF
% Builds aerodynamic coefficients and body-axis aero forces/moments
% from the current 6-DOF body state.
%
% state fields required:
%   u, v, w, p, q, r
%
% Outputs in aero:
%   Vt, alpha_deg, beta_deg
%   CXb, CYb, CZb
%   Cl, Cm, Cn
%   qbar

u = state.u;
v = state.v;
w = state.w;
p = state.p;
q = state.q;
r = state.r;

Vt = sqrt(u^2 + v^2 + w^2);
Veff = max(Vt, 1e-6);

alpha = atan2(w, u);                        % rad
beta  = asin(max(min(v / Veff, 1), -1));    % rad

alpha_deg = rad2deg(alpha);
beta_deg  = rad2deg(beta);

rho = state.rho;
qbar = 0.5 * rho * Vt^2;

% ----- Static coefficient lookup -----
BETA = db.aero.BETA;

CN0 = interp2(db.aero.ALFALCN1,   BETA, db.aero.LCN1,   alpha_deg, beta_deg, 'linear');
CA0 = interp2(db.aero.ALFALCA1,   BETA, db.aero.LCA1,   alpha_deg, beta_deg, 'linear');
CM0 = interp2(db.aero.ALFALCM1,   BETA, db.aero.LCM1,   alpha_deg, beta_deg, 'linear');
CY0 = interp2(db.aero.ALFALCY1,   BETA, db.aero.LCY1,   alpha_deg, beta_deg, 'linear');
Cl0 = interp2(db.aero.ALFALCSL1B, BETA, db.aero.LCSL1B, alpha_deg, beta_deg, 'linear');
Cn0 = interp2(db.aero.ALFALCSN1B, BETA, db.aero.LCSN1B, alpha_deg, beta_deg, 'linear');

% ----- Dynamic derivative lookup -----
p_hat = p * db.b / (2 * Veff);
q_hat = q * db.c / (2 * Veff);
r_hat = r * db.b / (2 * Veff);

CYp = interp1(db.aero.ALFADERIV,    db.aero.CYP, alpha_deg, 'linear', 'extrap');
CYr = interp1(db.aero.ALFADERIV,    db.aero.CYR, alpha_deg, 'linear', 'extrap');

Clp = interp1(db.aero.ALFADERIVCLP, db.aero.CLP, alpha_deg, 'linear', 'extrap');
Clr = interp1(db.aero.ALFADERIVCLR, db.aero.CLR, alpha_deg, 'linear', 'extrap');

Cmq = interp1(db.aero.ALFADERIVCMQ, db.aero.CMQ, alpha_deg, 'linear', 'extrap');

Cnp = interp1(db.aero.ALFADERIVCNP, db.aero.CNP, alpha_deg, 'linear', 'extrap');
Cnr = interp1(db.aero.ALFADERIVCNR, db.aero.CNR, alpha_deg, 'linear', 'extrap');

CY = CY0 + CYp * p_hat + CYr * r_hat;
Cl = Cl0 + Clp * p_hat + Clr * r_hat;
Cm = CM0 + Cmq * q_hat;
Cn = Cn0 + Cnp * p_hat + Cnr * r_hat;

% ----- Convert CA/CN to CL/CD -----
CL = CN0 * cos(alpha) - CA0 * sin(alpha);
CD = CN0 * sin(alpha) + CA0 * cos(alpha);

% ----- Wind-axis aero force coefficients -----
% x_w forward along velocity, y_w right, z_w down
CXw = -CD;
CYw = CY;
CZw = -CL;

% ----- Rotate wind-axis aero forces to body axes -----
R_w2b = dcm_wind_to_body(alpha, beta);
C_body = R_w2b * [CXw; CYw; CZw];

aero = struct();
aero.Vt        = Vt;
aero.alpha_deg = alpha_deg;
aero.beta_deg  = beta_deg;
aero.alpha_rad = alpha;
aero.beta_rad  = beta;
aero.qbar      = qbar;

aero.CL = CL;
aero.CD = CD;
aero.CY = CY;

aero.CXb = C_body(1);
aero.CYb = C_body(2);
aero.CZb = C_body(3);

aero.Cl = Cl;
aero.Cm = Cm;
aero.Cn = Cn;
end

function R = dcm_wind_to_body(alpha, beta)
% wind -> body using alpha, beta
ca = cos(alpha); sa = sin(alpha);
cb = cos(beta);  sb = sin(beta);

Rz_beta = [ cb -sb  0;
            sb  cb  0;
             0   0  1];

Ry_malpha = [ ca 0 -sa;
               0 1   0;
              sa 0  ca];

R = Ry_malpha * Rz_beta;
end