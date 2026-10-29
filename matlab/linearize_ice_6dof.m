function [A,B] = linearize_ice_6dof(x_trim, u_trim, db)
% LINEARIZE_ICE_6DOF
%
% Numerical linearisation of the 12-state 6-DOF model:
%
%   x = [u v w p q r phi theta psi xN yE zD]
%   xdot = ice_step_6dof(0, x, u, db)
%
% Inputs:
%   x_trim : 12x1 trim state vector
%   u_trim : struct with fields PLA, DELTAZ, DELTAW
%   db     : model database from ice_init()
%
% Outputs:
%   A : df/dx at trim
%   B : df/du at trim for inputs [PLA, DELTAZ, DELTAW]

% ---------------- checks ----------------
if ~isvector(x_trim) || numel(x_trim) ~= 12
    error('x_trim must be a 12x1 vector: [u v w p q r phi theta psi xN yE zD].');
end

req = {'PLA','DELTAZ','DELTAW'};
for k = 1:numel(req)
    if ~isfield(u_trim, req{k})
        error('u_trim must contain field "%s".', req{k});
    end
end

x_trim = x_trim(:);

nx = numel(x_trim);
A = zeros(nx,nx);

% ---------------- A matrix ----------------
for j = 1:nx
    x_plus  = x_trim;
    x_minus = x_trim;

    h = choose_step(x_trim(j));

    x_plus(j)  = x_plus(j)  + h;
    x_minus(j) = x_minus(j) - h;

    f_plus  = ice_step_6dof(0, x_plus,  u_trim, db);
    f_minus = ice_step_6dof(0, x_minus, u_trim, db);

    A(:,j) = (f_plus - f_minus) / (2*h);
end

% ---------------- B matrix ----------------
input_names = {'PLA','DELTAZ','DELTAW'};
nu = numel(input_names);
B = zeros(nx,nu);

for j = 1:nu
    name = input_names{j};

    u_plus  = u_trim;
    u_minus = u_trim;

    h = choose_step(u_trim.(name));

    u_plus.(name)  = u_plus.(name)  + h;
    u_minus.(name) = u_minus.(name) - h;

    f_plus  = ice_step_6dof(0, x_trim, u_plus,  db);
    f_minus = ice_step_6dof(0, x_trim, u_minus, db);

    B(:,j) = (f_plus - f_minus) / (2*h);
end

end

% =====================================================================
% step-size helper
% =====================================================================

function h = choose_step(val)
if ~isfinite(val)
    h = 1e-5;
else
    h = 1e-5 * max(1, abs(val));
end
end