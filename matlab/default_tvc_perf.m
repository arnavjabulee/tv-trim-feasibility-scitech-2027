function perf = default_tvc_perf()
% DEFAULT_TVC_PERF
%
% Static TVC performance terms
%   eta_tvc             scales turned thrust components / moment authority
%   k_loss              extra axial thrust loss with vector angle
%   delta_loss_ref_deg  normalisation angle used in the axial-loss model
%
% Dynamic TVC terms
%   delta_max_deg           hard gimbal angle limit [deg]
%   delta_dot_max_deg_s     gimbal slew-rate limit [deg/s]
%   tau_tvc_s               first-order TVC lag time constant [s]
%
% Dynamic throttle / thrust terms
%   tau_thrust_s            first-order lag on PLA state [s]
%   PLA_dot_max_per_s       PLA slew-rate limit [PLA units / s]

perf = struct();

% ---------------------------------
% Existing static performance terms
% ---------------------------------
perf.eta_tvc = 1.0;          % 1.0 = ideal TVC effectiveness
perf.k_loss  = 0.0;          % 0.0 = no extra axial thrust loss

% Keep this separate from the actuator hard stop.
% Your old model used delta_max_deg only as a loss normalisation angle.
perf.delta_loss_ref_deg = 25.0;   % [deg]

% ---------------------------------
% TVC actuator dynamics / limits
% ---------------------------------
perf.delta_max_deg = 3.0;          % hard gimbal angle limit [deg]
perf.delta_dot_max_deg_s = 5.0;    % gimbal slew-rate limit [deg/s]
perf.tau_tvc_s = 0.06;             % first-order actuator lag [s]

% ---------------------------------
% Throttle / thrust dynamics / limits
% ---------------------------------
perf.tau_thrust_s = 0.30;          % first-order lag on PLA [s]
perf.PLA_dot_max_per_s = inf;      % no PLA rate limit by default

end