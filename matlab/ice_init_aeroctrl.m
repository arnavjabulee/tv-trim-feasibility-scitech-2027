function db = ice_init_aeroctrl()
% ICE_INIT_AEROCTRL
%
% Initialises database for the aerodynamic-controls-only 6DOF branch.

load_ice_model_data();

db = struct();

% ============================================================
% Geometry / inertia
% ============================================================
db.S   = evalin('base','Sref');
db.b   = evalin('base','b');
db.c   = evalin('base','cbar');
db.m   = evalin('base','m');
db.g   = evalin('base','g');

db.xjFG = evalin('base','xjFG');
db.zjFG = evalin('base','zjFG');

db.xjDD = evalin('base','xjDD');
db.zjDD = evalin('base','zjDD');

db.Ixx = evalin('base','Ixx');
db.Iyy = evalin('base','Iyy');
db.Izz = evalin('base','Izz');
db.Ixz = evalin('base','Ixz');

% ============================================================
% Atmosphere
% ============================================================
db.atm.h_m = evalin('base','height(:)');
db.atm.rho = evalin('base','rho(:)');
db.atm.T   = evalin('base','temperature(:)');

% ============================================================
% Engine
% ============================================================
db.eng.XMACH  = evalin('base','XMACH(:)');
db.eng.ALT_ft = evalin('base','ALT(:)');
db.eng.PLA    = evalin('base','PLA(:)');
db.eng.FGtab  = evalin('base','FGtab');
db.eng.DDtab  = evalin('base','DDtab');

% ============================================================
% Baseline aero tables
% ============================================================
db.aero.BETA = evalin('base','BETA(:)''');

db.aero.LCN1   = evalin('base','LCN1');
db.aero.LCA1   = evalin('base','LCA1');
db.aero.LCM1   = evalin('base','LCM1');
db.aero.LCY1   = evalin('base','LCY1');
db.aero.LCSL1B = evalin('base','LCSL1B');
db.aero.LCSN1B = evalin('base','LCSN1B');

db.aero.ALFALCN1   = evalin('base','ALFA');
db.aero.ALFALCA1   = evalin('base','ALFALCA1');
db.aero.ALFALCM1   = evalin('base','ALFA');
db.aero.ALFALCY1   = evalin('base','ALFALCY1');
db.aero.ALFALCSL1B = evalin('base','ALFALCSL1B');
db.aero.ALFALCSN1B = evalin('base','ALFALCSN1B');

% ============================================================
% Dynamic derivatives
% ============================================================
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

% ============================================================
% Control-grid vectors
% ============================================================
db.aero.ctrl = struct();

db.aero.ctrl.ALFA       = evalin('base','ALFA(:)''');
db.aero.ctrl.ALFASHORT  = evalin('base','ALFASHORT(:)''');
db.aero.ctrl.ALFALCAD4  = evalin('base','ALFALCAD4(:)''');
db.aero.ctrl.ALFALCND4  = evalin('base','ALFALCND4(:)''');
db.aero.ctrl.BETA       = evalin('base','BETA(:)''');

db.aero.ctrl.IBLFL  = evalin('base','IBLFL(:)''');
db.aero.ctrl.OBLFL1 = evalin('base','OBLFL1(:)''');
db.aero.ctrl.OBLFL2 = evalin('base','OBLFL2(:)''');
db.aero.ctrl.ELEVL  = evalin('base','ELEVL(:)''');
db.aero.ctrl.SPLRL1 = evalin('base','SPLRL1(:)''');
db.aero.ctrl.SPLRL2 = evalin('base','SPLRL2(:)''');
db.aero.ctrl.SPLRL3 = evalin('base','SPLRL3(:)''');
db.aero.ctrl.SPLRL4 = evalin('base','SPLRL4(:)''');
db.aero.ctrl.PFLAP  = evalin('base','PFLAP(:)''');
db.aero.ctrl.CLAML  = evalin('base','CLAML(:)''');
db.aero.ctrl.AMTL   = evalin('base','AMTL(:)''');
db.aero.ctrl.DEFL   = evalin('base','DEFL(:)''');

db.aero.ctrl.IBLFR  = evalin('base','IBLFR(:)''');
db.aero.ctrl.OBLFR1 = evalin('base','OBLFR1(:)''');
db.aero.ctrl.OBLFR2 = evalin('base','OBLFR2(:)''');
db.aero.ctrl.ELEVR  = evalin('base','ELEVR(:)''');
db.aero.ctrl.SPLRR1 = evalin('base','SPLRR1(:)''');
db.aero.ctrl.SPLRR2 = evalin('base','SPLRR2(:)''');
db.aero.ctrl.SPLRR3 = evalin('base','SPLRR3(:)''');
db.aero.ctrl.SPLRR4 = evalin('base','SPLRR4(:)''');
db.aero.ctrl.CLAMR  = evalin('base','CLAMR(:)''');
db.aero.ctrl.AMTR   = evalin('base','AMTR(:)''');
db.aero.ctrl.DEFR   = evalin('base','DEFR(:)''');

% ============================================================
% Left-side increment tables
% ============================================================
db.aero.ctrl.LCAD1Ltab    = evalin('base','LCAD1Ltab');
db.aero.ctrl.LCND1Ltab    = evalin('base','LCND1Ltab');
db.aero.ctrl.LCMD1Ltab    = evalin('base','LCMD1Ltab');
db.aero.ctrl.LCYD1Ltab    = evalin('base','LCYD1Ltab');
db.aero.ctrl.LCSLD1LBtab  = evalin('base','LCSLD1LBtab');
db.aero.ctrl.LCSND1LBtab  = evalin('base','LCSND1LBtab');

db.aero.ctrl.LCAD2Ltab    = evalin('base','LCAD2Ltab');
db.aero.ctrl.LCND2Ltab    = evalin('base','LCND2Ltab');
db.aero.ctrl.LCMD2Ltab    = evalin('base','LCMD2Ltab');
db.aero.ctrl.LCYD2Ltab    = evalin('base','LCYD2Ltab');
db.aero.ctrl.LCSLD2LBtab  = evalin('base','LCSLD2LBtab');
db.aero.ctrl.LCSND2LBtab  = evalin('base','LCSND2LBtab');

db.aero.ctrl.LCAD3Ltab    = evalin('base','LCAD3Ltab');
db.aero.ctrl.LCND3Ltab    = evalin('base','LCND3Ltab');
db.aero.ctrl.LCMD3Ltab    = evalin('base','LCMD3Ltab');
db.aero.ctrl.LCYD3Ltab    = evalin('base','LCYD3Ltab');
db.aero.ctrl.LCSLD3LBtab  = evalin('base','LCSLD3LBtab');
db.aero.ctrl.LCSND3LBtab  = evalin('base','LCSND3LBtab');

db.aero.ctrl.LCAD4tab     = evalin('base','LCAD4tab');
db.aero.ctrl.LCND4tab     = evalin('base','LCND4tab');
db.aero.ctrl.LCMD4tab     = evalin('base','LCMD4tab');
db.aero.ctrl.LCYD4tab     = evalin('base','LCYD4tab');
db.aero.ctrl.LCSLD4Btab   = evalin('base','LCSLD4Btab');
db.aero.ctrl.LCSND4Btab   = evalin('base','LCSND4Btab');

db.aero.ctrl.LCAD5Ltab    = evalin('base','LCAD5Ltab');
db.aero.ctrl.LCND5Ltab    = evalin('base','LCND5Ltab');
db.aero.ctrl.LCMD5Ltab    = evalin('base','LCMD5Ltab');
db.aero.ctrl.LCYD5Ltab    = evalin('base','LCYD5Ltab');
db.aero.ctrl.LCSLD5LBtab  = evalin('base','LCSLD5LBtab');
db.aero.ctrl.LCSND5LBtab  = evalin('base','LCSND5LBtab');

db.aero.ctrl.LCAD6Ltab    = evalin('base','LCAD6Ltab');
db.aero.ctrl.LCND6Ltab    = evalin('base','LCND6Ltab');
db.aero.ctrl.LCMD6Ltab    = evalin('base','LCMD6Ltab');
db.aero.ctrl.LCYD6Ltab    = evalin('base','LCYD6Ltab');
db.aero.ctrl.LCSLD6LBtab  = evalin('base','LCSLD6LBtab');
db.aero.ctrl.LCSND6LBtab  = evalin('base','LCSND6LBtab');

db.aero.ctrl.LCAD8Ltab    = evalin('base','LCAD8Ltab');
db.aero.ctrl.LCND8Ltab    = evalin('base','LCND8Ltab');
db.aero.ctrl.LCMD8Ltab    = evalin('base','LCMD8Ltab');
db.aero.ctrl.LCYD8Ltab    = evalin('base','LCYD8Ltab');
db.aero.ctrl.LCSLD8LBtab  = evalin('base','LCSLD8LBtab');
db.aero.ctrl.LCSND8LBtab  = evalin('base','LCSND8LBtab');

db.aero.ctrl.LCAD9Ltab    = evalin('base','LCAD9Ltab');
db.aero.ctrl.LCND9Ltab    = evalin('base','LCND9Ltab');
db.aero.ctrl.LCMD9Ltab    = evalin('base','LCMD9Ltab');
db.aero.ctrl.LCYD9Ltab    = evalin('base','LCYD9Ltab');
db.aero.ctrl.LCSLD9LBtab  = evalin('base','LCSLD9LBtab');
db.aero.ctrl.LCSND9LBtab  = evalin('base','LCSND9LBtab');

% ============================================================
% Right-side increment tables
% ============================================================
db.aero.ctrl.LCAD11Rtab   = evalin('base','LCAD11Rtab');
db.aero.ctrl.LCND11Rtab   = evalin('base','LCND11Rtab');
db.aero.ctrl.LCMD11Rtab   = evalin('base','LCMD11Rtab');
db.aero.ctrl.LCYD11Rtab   = evalin('base','LCYD11Rtab');
db.aero.ctrl.LCSLD11RBtab = evalin('base','LCSLD11RBtab');
db.aero.ctrl.LCSND11RBtab = evalin('base','LCSND11RBtab');

db.aero.ctrl.LCAD12Rtab   = evalin('base','LCAD12Rtab');
db.aero.ctrl.LCND12Rtab   = evalin('base','LCND12Rtab');
db.aero.ctrl.LCMD12Rtab   = evalin('base','LCMD12Rtab');
db.aero.ctrl.LCYD12Rtab   = evalin('base','LCYD12Rtab');
db.aero.ctrl.LCSLD12RBtab = evalin('base','LCSLD12RBtab');
db.aero.ctrl.LCSND12RBtab = evalin('base','LCSND12RBtab');

db.aero.ctrl.LCAD13Rtab   = evalin('base','LCAD13Rtab');
db.aero.ctrl.LCND13Rtab   = evalin('base','LCND13Rtab');
db.aero.ctrl.LCMD13Rtab   = evalin('base','LCMD13Rtab');
db.aero.ctrl.LCYD13Rtab   = evalin('base','LCYD13Rtab');
db.aero.ctrl.LCSLD13RBtab = evalin('base','LCSLD13RBtab');
db.aero.ctrl.LCSND13RBtab = evalin('base','LCSND13RBtab');

db.aero.ctrl.LCAD15Rtab   = evalin('base','LCAD15Rtab');
db.aero.ctrl.LCND15Rtab   = evalin('base','LCND15Rtab');
db.aero.ctrl.LCMD15Rtab   = evalin('base','LCMD15Rtab');
db.aero.ctrl.LCYD15Rtab   = evalin('base','LCYD15Rtab');
db.aero.ctrl.LCSLD15RBtab = evalin('base','LCSLD15RBtab');
db.aero.ctrl.LCSND15RBtab = evalin('base','LCSND15RBtab');

db.aero.ctrl.LCAD16Rtab   = evalin('base','LCAD16Rtab');
db.aero.ctrl.LCND16Rtab   = evalin('base','LCND16Rtab');
db.aero.ctrl.LCMD16Rtab   = evalin('base','LCMD16Rtab');
db.aero.ctrl.LCYD16Rtab   = evalin('base','LCYD16Rtab');
db.aero.ctrl.LCSLD16RBtab = evalin('base','LCSLD16RBtab');
db.aero.ctrl.LCSND16RBtab = evalin('base','LCSND16RBtab');

db.aero.ctrl.LCAD18Rtab   = evalin('base','LCAD18Rtab');
db.aero.ctrl.LCND18Rtab   = evalin('base','LCND18Rtab');
db.aero.ctrl.CMD18Rtab    = evalin('base','LCMD18Rtab');
db.aero.ctrl.LCYD18Rtab   = evalin('base','LCYD18Rtab');
db.aero.ctrl.LCSLD18RBtab = evalin('base','LCSLD18RBtab');
db.aero.ctrl.LCSND18RBtab = evalin('base','LCSND18RBtab');

db.aero.ctrl.LCAD19Rtab   = evalin('base','LCAD19Rtab');
db.aero.ctrl.LCND19Rtab   = evalin('base','LCND19Rtab');
db.aero.ctrl.LCMD19Rtab   = evalin('base','LCMD19Rtab');
db.aero.ctrl.LCYD19Rtab   = evalin('base','LCYD19Rtab');
db.aero.ctrl.LCSLD19RBtab = evalin('base','LCSLD19RBtab');
db.aero.ctrl.LCSND19RBtab = evalin('base','LCSND19RBtab');

disp('ice_init_aeroctrl: OK')
end
