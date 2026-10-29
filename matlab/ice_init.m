function db = ice_init()

load_ice_model_data();

db = struct();

% Geometry
db.S   = evalin('base','Sref');
db.b   = evalin('base','b');
db.c   = evalin('base','cbar');
db.m   = evalin('base','m');
db.g   = evalin('base','g');

% --- Engine lever arms (from constants.m) ---
db.xjFG = evalin('base','xjFG');
db.zjFG = evalin('base','zjFG');

db.xjDD = evalin('base','xjDD');
db.zjDD = evalin('base','zjDD');


db.Ixx = evalin('base','Ixx');
db.Iyy = evalin('base','Iyy');
db.Izz = evalin('base','Izz');
db.Ixz = evalin('base','Ixz');

% Atmosphere
db.atm.h_m = evalin('base','height(:)');
db.atm.rho = evalin('base','rho(:)');
db.atm.T   = evalin('base','temperature(:)');

% Engine
db.eng.XMACH  = evalin('base','XMACH(:)');
db.eng.ALT_ft = evalin('base','ALT(:)');
db.eng.PLA    = evalin('base','PLA(:)');
db.eng.FGtab  = evalin('base','FGtab');
db.eng.DDtab  = evalin('base','DDtab');

% ---- Baseline aero tables ----
db.aero.BETA = evalin('base','BETA(:)''');

db.aero.LCN1   = evalin('base','LCN1');
db.aero.LCA1   = evalin('base','LCA1');
db.aero.LCM1   = evalin('base','LCM1');
db.aero.LCY1   = evalin('base','LCY1');
db.aero.LCSL1B = evalin('base','LCSL1B');
db.aero.LCSN1B = evalin('base','LCSN1B');

db.aero.ALFALCN1   = evalin('base','ALFA');        % 20
db.aero.ALFALCM1   = evalin('base','ALFA');        % 20
db.aero.ALFALCA1   = evalin('base','ALFALCA1');    % 26
db.aero.ALFALCY1   = evalin('base','ALFALCY1');    % 32
db.aero.ALFALCSL1B = evalin('base','ALFALCSL1B');  % 29
db.aero.ALFALCSN1B = evalin('base','ALFALCSN1B');  % 39

db.aero.ALFADERIV    = evalin('base','ALFADERIV');
db.aero.ALFADERIVCLP = evalin('base','ALFADERIVCLP');
db.aero.ALFADERIVCLR = evalin('base','ALFADERIVCLR');
db.aero.ALFADERIVCMQ = evalin('base','ALFADERIVCMQ');
db.aero.ALFADERIVCNP = evalin('base','ALFADERIVCNP');
db.aero.ALFADERIVCNR = evalin('base','ALFADERIVCNR');

db.aero.CYP = evalin('base','CYP');
db.aero.CYR = evalin('base','CYR');
db.aero.CLP = evalin('base','CLP');
db.aero.CLR = evalin('base','CLR');
db.aero.CMQ = evalin('base','CMQ');
db.aero.CNP = evalin('base','CNP');
db.aero.CNR = evalin('base','CNR');

disp('ice_init: OK')

end
