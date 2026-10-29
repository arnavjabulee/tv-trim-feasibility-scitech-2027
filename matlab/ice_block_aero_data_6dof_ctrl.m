function aero = ice_block_aero_data_6dof_ctrl(state, eff, db)
% ICE_BLOCK_AERO_DATA_6DOF_CTRL

u = state.u;
v = state.v;
w = state.w;
p = state.p;
q = state.q;
r = state.r;

Vt = sqrt(u^2 + v^2 + w^2);
Veff = max(Vt, 1e-6);

alpha = atan2(w, u);
beta  = asin(max(min(v / Veff, 1), -1));

alpha_deg = rad2deg(alpha);
beta_deg  = rad2deg(beta);

rho = state.rho;
qbar = 0.5 * rho * Vt^2;

% ------------------------------------------------------------
% Baseline static coefficient lookup
% ------------------------------------------------------------
BETA = db.aero.BETA;

CN0 = interpn(BETA, db.aero.ALFALCN1,   db.aero.LCN1,   beta_deg, alpha_deg, 'linear', 0);
CA0 = interpn(BETA, db.aero.ALFALCA1,   db.aero.LCA1,   beta_deg, alpha_deg, 'linear', 0);
CM0 = interpn(BETA, db.aero.ALFALCM1,   db.aero.LCM1,   beta_deg, alpha_deg, 'linear', 0);
CY0 = interpn(BETA, db.aero.ALFALCY1,   db.aero.LCY1,   beta_deg, alpha_deg, 'linear', 0);
Cl0 = interpn(BETA, db.aero.ALFALCSL1B, db.aero.LCSL1B, beta_deg, alpha_deg, 'linear', 0);
Cn0 = interpn(BETA, db.aero.ALFALCSN1B, db.aero.LCSN1B, beta_deg, alpha_deg, 'linear', 0);

% ------------------------------------------------------------
% Dynamic derivative lookup
% ------------------------------------------------------------
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

% ------------------------------------------------------------
% Read controls, defaulting to zero
% ------------------------------------------------------------
IBLFL = get_eff(eff, 'IBLFL', 0);
OBLFL = get_eff(eff, 'OBLFL', 0);
ELEVL = get_eff(eff, 'ELEVL', 0);
SPLRL = get_eff(eff, 'SPLRL', 0);
PFLAP = get_eff(eff, 'PFLAP', 0);
CLAML = get_eff(eff, 'CLAML', 0);
AMTL  = get_eff(eff, 'AMTL',  0);
DEFL  = get_eff(eff, 'DEFL',  0);

IBLFR = get_eff(eff, 'IBLFR', 0);
OBLFR = get_eff(eff, 'OBLFR', 0);
ELEVR = get_eff(eff, 'ELEVR', 0);
SPLRR = get_eff(eff, 'SPLRR', 0);
CLAMR = get_eff(eff, 'CLAMR', 0);
AMTR  = get_eff(eff, 'AMTR',  0);
DEFR  = get_eff(eff, 'DEFR',  0);

% ------------------------------------------------------------
% Control-surface coefficient increments
% ------------------------------------------------------------
dCA = 0;
dCN = 0;
dCM = 0;
dCY = 0;
dCl = 0;
dCn = 0;

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_1L(IBLFL, alpha_deg, beta_deg, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_2L(OBLFL, alpha_deg, IBLFL, beta_deg, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_3L(ELEVL, alpha_deg, SPLRL, OBLFL, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_4(PFLAP, alpha_deg, SPLRL, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_5L(beta_deg, alpha_deg, SPLRL, CLAML, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_6L(SPLRL, alpha_deg, OBLFL, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_8L(AMTL, alpha_deg, beta_deg, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_9L(DEFL, alpha_deg, SPLRL, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_11R(IBLFR, alpha_deg, beta_deg, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_12R(OBLFR, alpha_deg, IBLFR, beta_deg, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_13R(ELEVR, alpha_deg, SPLRR, OBLFR, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_15R(beta_deg, alpha_deg, SPLRR, CLAMR, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_16R(SPLRR, alpha_deg, OBLFR, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_18R(AMTR, alpha_deg, beta_deg, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

[dca,dcn,dcm,dcy,dcl,dcnyaw] = inc_block_19R(DEFR, alpha_deg, SPLRR, db);
[dCA,dCN,dCM,dCY,dCl,dCn] = add_inc(dCA,dCN,dCM,dCY,dCl,dCn,dca,dcn,dcm,dcy,dcl,dcnyaw);

% ------------------------------------------------------------
% Combine baseline + increments + dynamic derivatives
% ------------------------------------------------------------
CA = CA0 + dCA;
CN = CN0 + dCN;
CM = CM0 + dCM;

CY = CY0 + dCY + CYp * p_hat + CYr * r_hat;
Cl = Cl0 + dCl + Clp * p_hat + Clr * r_hat;
Cm = CM  + Cmq * q_hat;
Cn = Cn0 + dCn + Cnp * p_hat + Cnr * r_hat;

% ------------------------------------------------------------
% Convert CA/CN to CL/CD
% ------------------------------------------------------------
CL = CN * cos(alpha) - CA * sin(alpha);
CD = CN * sin(alpha) + CA * cos(alpha);

CXw = -CD;
CYw = CY;
CZw = -CL;

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

aero.debug = struct( ...
    'dCA', dCA, ...
    'dCN', dCN, ...
    'dCM', dCM, ...
    'dCY', dCY, ...
    'dCl', dCl, ...
    'dCn', dCn);
end

% ============================================================
% Helpers
% ============================================================

function [CA,CN,CM,CY,Cl,Cn] = add_inc(CA,CN,CM,CY,Cl,Cn,dCA,dCN,dCM,dCY,dCl,dCn)
CA = CA + dCA;
CN = CN + dCN;
CM = CM + dCM;
CY = CY + dCY;
Cl = Cl + dCl;
Cn = Cn + dCn;
end

function val = get_eff(eff, name, defaultVal)
if isfield(eff, name) && ~isempty(eff.(name)) && isfinite(eff.(name))
    val = eff.(name);
else
    val = defaultVal;
end
end

function R = dcm_wind_to_body(alpha, beta)
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

% ============================================================
% Increment blocks
% ============================================================

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_1L(iblfl, alpha, beta, db)
dCA = interpn(db.aero.ctrl.IBLFL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCAD1Ltab, iblfl, alpha, beta, 'linear', 0);
dCN = interpn(db.aero.ctrl.IBLFL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCND1Ltab, iblfl, alpha, beta, 'linear', 0);
dCM = interpn(db.aero.ctrl.IBLFL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCMD1Ltab, iblfl, alpha, beta, 'linear', 0);
dCY = interpn(db.aero.ctrl.IBLFL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCYD1Ltab, iblfl, alpha, beta, 'linear', 0);
dCl = interpn(db.aero.ctrl.IBLFL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSLD1LBtab, iblfl, alpha, beta, 'linear', 0);
dCn = interpn(db.aero.ctrl.IBLFL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSND1LBtab, iblfl, alpha, beta, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_2L(oblfl, alpha, iblfl, beta, db)
dCA = interpn(db.aero.ctrl.OBLFL1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFL, db.aero.ctrl.BETA, db.aero.ctrl.LCAD2Ltab, oblfl, alpha, iblfl, beta, 'linear', 0);
dCN = interpn(db.aero.ctrl.OBLFL1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFL, db.aero.ctrl.BETA, db.aero.ctrl.LCND2Ltab, oblfl, alpha, iblfl, beta, 'linear', 0);
dCM = interpn(db.aero.ctrl.OBLFL1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFL, db.aero.ctrl.BETA, db.aero.ctrl.LCMD2Ltab, oblfl, alpha, iblfl, beta, 'linear', 0);
dCY = interpn(db.aero.ctrl.OBLFL1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFL, db.aero.ctrl.BETA, db.aero.ctrl.LCYD2Ltab, oblfl, alpha, iblfl, beta, 'linear', 0);
dCl = interpn(db.aero.ctrl.OBLFL1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFL, db.aero.ctrl.BETA, db.aero.ctrl.LCSLD2LBtab, oblfl, alpha, iblfl, beta, 'linear', 0);
dCn = interpn(db.aero.ctrl.OBLFL1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFL, db.aero.ctrl.BETA, db.aero.ctrl.LCSND2LBtab, oblfl, alpha, iblfl, beta, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_3L(elevl, alpha, splrl, oblfl, db)
dCA = interpn(db.aero.ctrl.ELEVL, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL1, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCAD3Ltab, elevl, alpha, splrl, oblfl, 'linear', 0);
dCN = interpn(db.aero.ctrl.ELEVL, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL1, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCND3Ltab, elevl, alpha, splrl, oblfl, 'linear', 0);
dCM = interpn(db.aero.ctrl.ELEVL, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL1, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCMD3Ltab, elevl, alpha, splrl, oblfl, 'linear', 0);
dCY = interpn(db.aero.ctrl.ELEVL, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL1, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCYD3Ltab, elevl, alpha, splrl, oblfl, 'linear', 0);
dCl = interpn(db.aero.ctrl.ELEVL, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL1, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCSLD3LBtab, elevl, alpha, splrl, oblfl, 'linear', 0);
dCn = interpn(db.aero.ctrl.ELEVL, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL1, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCSND3LBtab, elevl, alpha, splrl, oblfl, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_4(pflap, alpha, splrl, db)
dCA = interpn(db.aero.ctrl.PFLAP, db.aero.ctrl.ALFALCAD4, db.aero.ctrl.SPLRL2, db.aero.ctrl.LCAD4tab, pflap, alpha, splrl, 'linear', 0);
dCN = interpn(db.aero.ctrl.PFLAP, db.aero.ctrl.ALFALCND4, db.aero.ctrl.SPLRL2, db.aero.ctrl.LCND4tab, pflap, alpha, splrl, 'linear', 0);
dCM = interpn(db.aero.ctrl.PFLAP, db.aero.ctrl.ALFA,       db.aero.ctrl.SPLRL2, db.aero.ctrl.LCMD4tab, pflap, alpha, splrl, 'linear', 0);
dCY = interpn(db.aero.ctrl.PFLAP, db.aero.ctrl.ALFA,       db.aero.ctrl.SPLRL2, db.aero.ctrl.LCYD4tab, pflap, alpha, splrl, 'linear', 0);
dCl = interpn(db.aero.ctrl.PFLAP, db.aero.ctrl.ALFA,       db.aero.ctrl.SPLRL2, db.aero.ctrl.LCSLD4Btab, pflap, alpha, splrl, 'linear', 0);
dCn = interpn(db.aero.ctrl.PFLAP, db.aero.ctrl.ALFA,       db.aero.ctrl.SPLRL2, db.aero.ctrl.LCSND4Btab, pflap, alpha, splrl, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_5L(beta, alpha, splrl, claml, db)
dCA = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL2, db.aero.ctrl.CLAML, db.aero.ctrl.LCAD5Ltab, beta, alpha, splrl, claml, 'linear', 0);
dCN = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL2, db.aero.ctrl.CLAML, db.aero.ctrl.LCND5Ltab, beta, alpha, splrl, claml, 'linear', 0);
dCM = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL2, db.aero.ctrl.CLAML, db.aero.ctrl.LCMD5Ltab, beta, alpha, splrl, claml, 'linear', 0);
dCY = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL2, db.aero.ctrl.CLAML, db.aero.ctrl.LCYD5Ltab, beta, alpha, splrl, claml, 'linear', 0);
dCl = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL2, db.aero.ctrl.CLAML, db.aero.ctrl.LCSLD5LBtab, beta, alpha, splrl, claml, 'linear', 0);
dCn = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRL2, db.aero.ctrl.CLAML, db.aero.ctrl.LCSND5LBtab, beta, alpha, splrl, claml, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_6L(splrl, alpha, oblfl, db)
dCA = interpn(db.aero.ctrl.SPLRL3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCAD6Ltab, splrl, alpha, oblfl, 'linear', 0);
dCN = interpn(db.aero.ctrl.SPLRL3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCND6Ltab, splrl, alpha, oblfl, 'linear', 0);
dCM = interpn(db.aero.ctrl.SPLRL3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCMD6Ltab, splrl, alpha, oblfl, 'linear', 0);
dCY = interpn(db.aero.ctrl.SPLRL3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCYD6Ltab, splrl, alpha, oblfl, 'linear', 0);
dCl = interpn(db.aero.ctrl.SPLRL3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCSLD6LBtab, splrl, alpha, oblfl, 'linear', 0);
dCn = interpn(db.aero.ctrl.SPLRL3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFL2, db.aero.ctrl.LCSND6LBtab, splrl, alpha, oblfl, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_8L(amtl, alpha, beta, db)
dCA = interpn(db.aero.ctrl.AMTL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCAD8Ltab, amtl, alpha, beta, 'linear', 0);
dCN = interpn(db.aero.ctrl.AMTL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCND8Ltab, amtl, alpha, beta, 'linear', 0);
dCM = interpn(db.aero.ctrl.AMTL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCMD8Ltab, amtl, alpha, beta, 'linear', 0);
dCY = interpn(db.aero.ctrl.AMTL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCYD8Ltab, amtl, alpha, beta, 'linear', 0);
dCl = interpn(db.aero.ctrl.AMTL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSLD8LBtab, amtl, alpha, beta, 'linear', 0);
dCn = interpn(db.aero.ctrl.AMTL, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSND8LBtab, amtl, alpha, beta, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_9L(defl, alpha, splrl, db)
dCA = interpn(db.aero.ctrl.DEFL, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRL4, db.aero.ctrl.LCAD9Ltab, defl, alpha, splrl, 'linear', 0);
dCN = interpn(db.aero.ctrl.DEFL, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRL4, db.aero.ctrl.LCND9Ltab, defl, alpha, splrl, 'linear', 0);
dCM = interpn(db.aero.ctrl.DEFL, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRL4, db.aero.ctrl.LCMD9Ltab, defl, alpha, splrl, 'linear', 0);
dCY = interpn(db.aero.ctrl.DEFL, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRL4, db.aero.ctrl.LCYD9Ltab, defl, alpha, splrl, 'linear', 0);
dCl = interpn(db.aero.ctrl.DEFL, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRL4, db.aero.ctrl.LCSLD9LBtab, defl, alpha, splrl, 'linear', 0);
dCn = interpn(db.aero.ctrl.DEFL, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRL4, db.aero.ctrl.LCSND9LBtab, defl, alpha, splrl, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_11R(iblfr, alpha, beta, db)
dCA = interpn(db.aero.ctrl.IBLFR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCAD11Rtab, iblfr, alpha, beta, 'linear', 0);
dCN = interpn(db.aero.ctrl.IBLFR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCND11Rtab, iblfr, alpha, beta, 'linear', 0);
dCM = interpn(db.aero.ctrl.IBLFR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCMD11Rtab, iblfr, alpha, beta, 'linear', 0);
dCY = interpn(db.aero.ctrl.IBLFR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCYD11Rtab, iblfr, alpha, beta, 'linear', 0);
dCl = interpn(db.aero.ctrl.IBLFR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSLD11RBtab, iblfr, alpha, beta, 'linear', 0);
dCn = interpn(db.aero.ctrl.IBLFR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSND11RBtab, iblfr, alpha, beta, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_12R(oblfr, alpha, iblfr, beta, db)
dCA = interpn(db.aero.ctrl.OBLFR1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFR, db.aero.ctrl.BETA, db.aero.ctrl.LCAD12Rtab, oblfr, alpha, iblfr, beta, 'linear', 0);
dCN = interpn(db.aero.ctrl.OBLFR1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFR, db.aero.ctrl.BETA, db.aero.ctrl.LCND12Rtab, oblfr, alpha, iblfr, beta, 'linear', 0);
dCM = interpn(db.aero.ctrl.OBLFR1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFR, db.aero.ctrl.BETA, db.aero.ctrl.LCMD12Rtab, oblfr, alpha, iblfr, beta, 'linear', 0);
dCY = interpn(db.aero.ctrl.OBLFR1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFR, db.aero.ctrl.BETA, db.aero.ctrl.LCYD12Rtab, oblfr, alpha, iblfr, beta, 'linear', 0);
dCl = interpn(db.aero.ctrl.OBLFR1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFR, db.aero.ctrl.BETA, db.aero.ctrl.LCSLD12RBtab, oblfr, alpha, iblfr, beta, 'linear', 0);
dCn = interpn(db.aero.ctrl.OBLFR1, db.aero.ctrl.ALFA, db.aero.ctrl.IBLFR, db.aero.ctrl.BETA, db.aero.ctrl.LCSND12RBtab, oblfr, alpha, iblfr, beta, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_13R(elevr, alpha, splrr, oblfr, db)
dCA = interpn(db.aero.ctrl.ELEVR, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR1, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCAD13Rtab, elevr, alpha, splrr, oblfr, 'linear', 0);
dCN = interpn(db.aero.ctrl.ELEVR, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR1, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCND13Rtab, elevr, alpha, splrr, oblfr, 'linear', 0);
dCM = interpn(db.aero.ctrl.ELEVR, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR1, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCMD13Rtab, elevr, alpha, splrr, oblfr, 'linear', 0);
dCY = interpn(db.aero.ctrl.ELEVR, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR1, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCYD13Rtab, elevr, alpha, splrr, oblfr, 'linear', 0);
dCl = interpn(db.aero.ctrl.ELEVR, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR1, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCSLD13RBtab, elevr, alpha, splrr, oblfr, 'linear', 0);
dCn = interpn(db.aero.ctrl.ELEVR, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR1, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCSND13RBtab, elevr, alpha, splrr, oblfr, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_15R(beta, alpha, splrr, clamr, db)
dCA = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR2, db.aero.ctrl.CLAMR, db.aero.ctrl.LCAD15Rtab, beta, alpha, splrr, clamr, 'linear', 0);
dCN = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR2, db.aero.ctrl.CLAMR, db.aero.ctrl.LCND15Rtab, beta, alpha, splrr, clamr, 'linear', 0);
dCM = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR2, db.aero.ctrl.CLAMR, db.aero.ctrl.LCMD15Rtab, beta, alpha, splrr, clamr, 'linear', 0);
dCY = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR2, db.aero.ctrl.CLAMR, db.aero.ctrl.LCYD15Rtab, beta, alpha, splrr, clamr, 'linear', 0);
dCl = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR2, db.aero.ctrl.CLAMR, db.aero.ctrl.LCSLD15RBtab, beta, alpha, splrr, clamr, 'linear', 0);
dCn = interpn(db.aero.ctrl.BETA, db.aero.ctrl.ALFA, db.aero.ctrl.SPLRR2, db.aero.ctrl.CLAMR, db.aero.ctrl.LCSND15RBtab, beta, alpha, splrr, clamr, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_16R(splrr, alpha, oblfr, db)
dCA = interpn(db.aero.ctrl.SPLRR3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCAD16Rtab, splrr, alpha, oblfr, 'linear', 0);
dCN = interpn(db.aero.ctrl.SPLRR3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCND16Rtab, splrr, alpha, oblfr, 'linear', 0);
dCM = interpn(db.aero.ctrl.SPLRR3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCMD16Rtab, splrr, alpha, oblfr, 'linear', 0);
dCY = interpn(db.aero.ctrl.SPLRR3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCYD16Rtab, splrr, alpha, oblfr, 'linear', 0);
dCl = interpn(db.aero.ctrl.SPLRR3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCSLD16RBtab, splrr, alpha, oblfr, 'linear', 0);
dCn = interpn(db.aero.ctrl.SPLRR3, db.aero.ctrl.ALFA, db.aero.ctrl.OBLFR2, db.aero.ctrl.LCSND16RBtab, splrr, alpha, oblfr, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_18R(amtr, alpha, beta, db)
dCA = interpn(db.aero.ctrl.AMTR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCAD18Rtab, amtr, alpha, beta, 'linear', 0);
dCN = interpn(db.aero.ctrl.AMTR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCND18Rtab, amtr, alpha, beta, 'linear', 0);
dCM = interpn(db.aero.ctrl.AMTR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.CMD18Rtab, amtr, alpha, beta, 'linear', 0);
dCY = interpn(db.aero.ctrl.AMTR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCYD18Rtab, amtr, alpha, beta, 'linear', 0);
dCl = interpn(db.aero.ctrl.AMTR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSLD18RBtab, amtr, alpha, beta, 'linear', 0);
dCn = interpn(db.aero.ctrl.AMTR, db.aero.ctrl.ALFA, db.aero.ctrl.BETA, db.aero.ctrl.LCSND18RBtab, amtr, alpha, beta, 'linear', 0);
end

function [dCA,dCN,dCM,dCY,dCl,dCn] = inc_block_19R(defr, alpha, splrr, db)
dCA = interpn(db.aero.ctrl.DEFR, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRR4, db.aero.ctrl.LCAD19Rtab, defr, alpha, splrr, 'linear', 0);
dCN = interpn(db.aero.ctrl.DEFR, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRR4, db.aero.ctrl.LCND19Rtab, defr, alpha, splrr, 'linear', 0);
dCM = interpn(db.aero.ctrl.DEFR, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRR4, db.aero.ctrl.LCMD19Rtab, defr, alpha, splrr, 'linear', 0);
dCY = interpn(db.aero.ctrl.DEFR, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRR4, db.aero.ctrl.LCYD19Rtab, defr, alpha, splrr, 'linear', 0);
dCl = interpn(db.aero.ctrl.DEFR, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRR4, db.aero.ctrl.LCSLD19RBtab, defr, alpha, splrr, 'linear', 0);
dCn = interpn(db.aero.ctrl.DEFR, db.aero.ctrl.ALFASHORT, db.aero.ctrl.SPLRR4, db.aero.ctrl.LCSND19RBtab, defr, alpha, splrr, 'linear', 0);
end