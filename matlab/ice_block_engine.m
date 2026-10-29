% function eng = ice_block_engine(z, Vt, PLA_cmd, db)
% 
% h_m = max(0,-z);
% 
% rho = interp1(db.atm.h_m, db.atm.rho, h_m, 'linear','extrap');
% rho = max(rho,0);
% 
% T = interp1(db.atm.h_m, db.atm.T, h_m, 'linear','extrap');
% 
% gamma = 1.4;
% R = 287.05;
% a = sqrt(gamma*R*T);
% mach = max(0, Vt/max(a,1e-6));
% 
% h_ft = h_m*3.28084;
% 
% mach_q = min(max(mach, min(db.eng.XMACH)), max(db.eng.XMACH));
% alt_q  = min(max(h_ft, min(db.eng.ALT_ft)), max(db.eng.ALT_ft));
% pla_q  = min(max(PLA_cmd, min(db.eng.PLA)), max(db.eng.PLA));
% 
% Fg = griddedInterpolant({db.eng.ALT_ft, db.eng.XMACH, db.eng.PLA}, db.eng.FGtab,'linear','nearest');
% FG = Fg(alt_q, mach_q, pla_q);
% 
% Fd = griddedInterpolant({db.eng.ALT_ft, db.eng.XMACH, db.eng.PLA}, db.eng.DDtab,'linear','nearest');
% DD = Fd(alt_q, mach_q, pla_q);
% 
% eng = struct('rho',rho,'FG',FG,'DD',DD);
% end

function eng = ice_block_engine(z, Vt, PLA_cmd, db)

h_m = max(0, -z);

rho = interp1(db.atm.h_m, db.atm.rho, h_m, 'linear', 'extrap');
rho = max(rho, 0);

T = interp1(db.atm.h_m, db.atm.T, h_m, 'linear', 'extrap');

gamma = 1.4;
R = 287.05;
a = sqrt(gamma * R * T);
mach = max(0, Vt / max(a, 1e-6));

h_ft = h_m * 3.28084;

% Keep altitude and Mach inside known table range
mach_q = min(max(mach, min(db.eng.XMACH)), max(db.eng.XMACH));
alt_q  = min(max(h_ft, min(db.eng.ALT_ft)), max(db.eng.ALT_ft));

pla_grid = db.eng.PLA(:);

FG = local_interp_extrap_pla(db.eng.FGtab, db, alt_q, mach_q, PLA_cmd, pla_grid);
DD = local_interp_extrap_pla(db.eng.DDtab, db, alt_q, mach_q, PLA_cmd, pla_grid);

% Optional physical floor
FG = max(FG, 0);
DD = max(DD, 0);

eng = struct('rho', rho, 'FG', FG, 'DD', DD);
end

% =====================================================================
% LOCAL HELPERS
% =====================================================================

function val = local_interp_extrap_pla(tab, db, alt_q, mach_q, pla_cmd, pla_grid)

% Inside table: normal interpolation
if pla_cmd >= pla_grid(1) && pla_cmd <= pla_grid(end)
    F = griddedInterpolant({db.eng.ALT_ft, db.eng.XMACH, db.eng.PLA}, ...
        tab, 'linear', 'nearest');
    val = F(alt_q, mach_q, pla_cmd);
    return
end

% Below table: linear extrapolation from first two PLA points
if pla_cmd < pla_grid(1)
    p1 = pla_grid(1);
    p2 = pla_grid(2);

    v1 = interpn(db.eng.ALT_ft, db.eng.XMACH, db.eng.PLA, ...
        tab, alt_q, mach_q, p1, 'linear');
    v2 = interpn(db.eng.ALT_ft, db.eng.XMACH, db.eng.PLA, ...
        tab, alt_q, mach_q, p2, 'linear');

    slope = (v2 - v1) / (p2 - p1);
    val = v1 + slope * (pla_cmd - p1);
    return
end

% Above table: linear extrapolation from last two PLA points
p1 = pla_grid(end-1);
p2 = pla_grid(end);

v1 = interpn(db.eng.ALT_ft, db.eng.XMACH, db.eng.PLA, ...
    tab, alt_q, mach_q, p1, 'linear');
v2 = interpn(db.eng.ALT_ft, db.eng.XMACH, db.eng.PLA, ...
    tab, alt_q, mach_q, p2, 'linear');

slope = (v2 - v1) / (p2 - p1);
val = v2 + slope * (pla_cmd - p2);

end