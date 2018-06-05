function [F, M] = controller(t, state, des_state, params)
%CONTROLLER  Controller for the quadrotor
%
%   state: The current state of the robot with the following fields:
%   state.pos = [x; y; z], state.vel = [x_dot; y_dot; z_dot],
%   state.rot = [phi; theta; psi], state.omega = [p; q; r]
%
%   des_state: The desired states are:
%   des_state.pos = [x; y; z], des_state.vel = [x_dot; y_dot; z_dot],
%   des_state.acc = [x_ddot; y_ddot; z_ddot], des_state.yaw,
%   des_state.yawdot
%
%   params: robot parameters


error = des_state.pos - state.pos;
% if norm(des_state.vel)
%     % tangent vector
%     t_unit = des_state.vel / norm(des_state.vel);
%     % normal vector
%     n      = des_state.acc / norm(des_state.vel) - dot(des_state.vel,des_state.acc)*des_state.vel/power(norm(des_state.vel),1.5);
%     n_unit = zeros(3,1);
%     if norm(n)
%         n_unit = n / norm(n);
%     end
%     %binormal vector
%     b      = cross(t_unit, n);
%     b_unit = zeros(3,1);
%     if norm(b)
%         b_unit = b / norm(b);
%     end
%     error = dot(error, n_unit)*n_unit + dot(error, b_unit)*b_unit;
% end    

omega_z = 50;
kv_z = 1.4*omega_z;
kp_z = omega_z*omega_z;
F = params.mass*(params.gravity + des_state.acc(3)  + kv_z*(des_state.vel(3)-state.vel(3)) + kp_z*error(3));

omega_x = 5;
kv_x = 1.4*omega_x;
kp_x = omega_x*omega_x;
r_des1 = des_state.acc(1)  + kv_x*(des_state.vel(1)-state.vel(1)) + kp_x*error(1);

omega_y = 5;
kv_y = 1.4*omega_y;
kp_y = omega_y*omega_y;
r_des2 = des_state.acc(2)  + kv_y*(des_state.vel(2)-state.vel(2)) + kp_y*error(2);

phi_des   = 1/params.gravity*(r_des1*sin(state.rot(3)) - r_des2*cos(state.rot(3)));
theta_des = 1/params.gravity*(r_des1*cos(state.rot(3)) + r_des2*sin(state.rot(3)));

% Moment
omega_phi   = 50;
omega_theta = 50;
omega_psi   = 50;

kv_phi   = 1.4*omega_phi  ;
kv_theta = 1.4*omega_theta;
kv_psi   = 1.4*omega_psi;

kp_phi   = omega_phi  *omega_phi  ;
kp_theta = omega_theta*omega_theta;
kp_psi   = omega_psi  *omega_psi  ;

ref1  = kv_phi  *                 -state.omega(1) + kp_phi  *(phi_des         -state.rot(1));
ref2  = kv_theta*                 -state.omega(2) + kp_theta*(theta_des       -state.rot(2));
ref3  = kv_psi  *(des_state.yawdot-state.omega(3))+ kp_psi  *(des_state.yawdot-state.rot(3));

omega_askew = [ 0             , -state.omega(3),  state.omega(2)
                state.omega(3),  0             , -state.omega(1)
               -state.omega(2),  state.omega(1),  0             ];
M(:,1) = params.I*[ref1; ref2; ref3] + omega_askew*params.I*state.omega;
%M = zeros(3,1);
% =================== Your code ends here ===================

end
