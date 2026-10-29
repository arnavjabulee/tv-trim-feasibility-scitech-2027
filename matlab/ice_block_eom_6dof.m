function xdot = ice_block_eom_6dof(state, FM, db)
% ICE_BLOCK_EOM_6DOF
% Full rigid-body 6-DOF equations using body velocities and Euler angles.
%
% State order:
%   x = [u v w p q r phi theta psi xN yE zD]

u = state.u;
v = state.v;
w = state.w;

p = state.p;
q = state.q;
r = state.r;

phi   = state.phi;
theta = state.theta;
psi   = state.psi;

Fb = FM.Fb;
Mb = FM.Mb;

m = db.m;
g = db.g;

% ----- Translational dynamics -----
uvw = [u; v; w];
omega = [p; q; r];

% Gravity in body axes (NED frame, z down)
g_b = [-g * sin(theta);
        g * sin(phi) * cos(theta);
        g * cos(phi) * cos(theta)];

uvw_dot = Fb / m - cross(omega, uvw) + g_b;

% ----- Rotational dynamics -----
I = [ db.Ixx,    0,    -db.Ixz;
         0,   db.Iyy,    0;
      -db.Ixz,   0,   db.Izz ];

omega_dot = I \ (Mb - cross(omega, I * omega));

% ----- Euler angle kinematics -----
cphi = cos(phi);   sphi = sin(phi);
cth  = cos(theta); sth  = sin(theta);
tth  = tan(theta);

E = [1, sphi*tth,  cphi*tth;
     0, cphi,     -sphi;
     0, sphi/cth,  cphi/cth];

euler_dot = E * omega;

% ----- Position kinematics (body -> NED) -----
cpsi = cos(psi); spsi = sin(psi);

R_b2n = [ cth*cpsi,  sphi*sth*cpsi - cphi*spsi,  cphi*sth*cpsi + sphi*spsi;
          cth*spsi,  sphi*sth*spsi + cphi*cpsi,  cphi*sth*spsi - sphi*cpsi;
             -sth,                sphi*cth,                 cphi*cth      ];

pos_dot = R_b2n * uvw;

xdot = [uvw_dot; omega_dot; euler_dot; pos_dot];
end