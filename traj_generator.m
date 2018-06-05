function [ desired_state ] = traj_generator(t, state, waypoints)
% TRAJ_GENERATOR: Generate the trajectory passing through all
% positions listed in the waypoints list
%
% t,state: time and current state (same variable as "state" in controller)
% that you may use for computing desired_state
%
% waypoints: The 3xP matrix listing all the points you much visited in order
% along the generated trajectory
%
% desired_state: Contains all the information that is passed to the
% controller for generating inputs for the quadrotor
%
persistent waypoints0 traj_time d0 alpha_x alpha_y alpha_z
if nargin > 2
    d = waypoints(:,2:end) - waypoints(:,1:end-1);
    d0 = 1* sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);
    traj_time = [0, cumsum(d0)];
    waypoints0 = waypoints;
    
    n = size(waypoints, 2) -1;
    
    A = zeros(8*n);
    % location constraints
    A(  1  :  n  ,  1:n  ) = eye(n);
    A(n+1  :2*n  ,   :   ) = repmat(eye(n),1,8);
    % speed constraints
    temp = eye(n) - diag(ones(1,n-1), 1);
    A(2*n+3:3*n+1,  n+1:end) = [prod(1:1)*temp(1:n-1,:), prod(2:2)*eye(n-1,n), prod(3:3)*eye(n-1,n), prod(4:4)*eye(n-1,n), prod(5:5)*eye(n-1,n), prod(6:6)*eye(n-1,n), prod(7:7)*eye(n-1,n)];
    % acceleration constraints
    A(3*n+2:4*n  ,2*n+1:end) = [prod(1:2)*temp(1:n-1,:), prod(2:3)*eye(n-1,n), prod(3:4)*eye(n-1,n), prod(4:5)*eye(n-1,n), prod(5:6)*eye(n-1,n), prod(6:7)*eye(n-1,n)];
    % Jerk constraints
    A(4*n+1:5*n-1,3*n+1:end) = [prod(1:3)*temp(1:n-1,:), prod(2:4)*eye(n-1,n), prod(3:5)*eye(n-1,n), prod(4:6)*eye(n-1,n), prod(5:7)*eye(n-1,n)];
    % snap constraints
    A(5*n  :6*n-2,4*n+1:end) = [prod(1:4)*temp(1:n-1,:), prod(2:5)*eye(n-1,n), prod(3:6)*eye(n-1,n), prod(4:7)*eye(n-1,n)];
    % dsnap constraints
    A(6*n-1:7*n-3,5*n+1:end) = [prod(1:5)*temp(1:n-1,:), prod(2:6)*eye(n-1,n), prod(3:7)*eye(n-1,n)];
    % ddsnap constraints
    A(7*n-2:8*n-4,6*n+1:end) = [prod(1:6)*temp(1:n-1,:), prod(2:7)*eye(n-1,n)];
        
    temp = zeros(1,n);
    temp(n) = 1;
    % initial and end speed constraints
    A(2*n+1      ,  n+1    ) = 1;
    A(2*n+2      ,  n+1:end) = [prod(1:1)*temp, prod(2:2)*temp, prod(3:3)*temp, prod(4:4)*temp, prod(5:5)*temp, prod(6:6)*temp, prod(7:7)*temp];
    % initial and end acceleration constraints
    A(8*n-3      ,2*n+1) = 1;
    A(8*n-2      ,2*n+1:end) = [prod(1:2)*temp, prod(2:3)*temp, prod(3:4)*temp, prod(4:5)*temp, prod(5:6)*temp, prod(6:7)*temp];
    % initial and end snap constraints
    A(8*n-1      ,3*n+1) = 1;
    A(8*n        ,3*n+1:end) = [prod(1:3)*temp, prod(2:4)*temp, prod(3:5)*temp, prod(4:6)*temp, prod(5:7)*temp];
    
    b_x = zeros(8*n,1);
    b_y = zeros(8*n,1);
    b_z = zeros(8*n,1);
    
    b_x(1  :  n) = waypoints(1,1:n  );
    b_x(n+1:2*n) = waypoints(1,2:n+1);
    
    b_y(1  :  n) = waypoints(2,1:n  );
    b_y(n+1:2*n) = waypoints(2,2:n+1);
    
    b_z(1  :  n) = waypoints(3,1:n  );
    b_z(n+1:2*n) = waypoints(3,2:n+1);
    
    alpha_x = A\b_x;
    alpha_y = A\b_y;
    alpha_z = A\b_z;
    
    alpha_x = fliplr(reshape(alpha_x,n,[]));
    alpha_y = fliplr(reshape(alpha_y,n,[]));
    alpha_z = fliplr(reshape(alpha_z,n,[]));
else
    if(t > traj_time(end))
        desired_state.pos = waypoints0(:,end);
        desired_state.vel = zeros(3,1);
        desired_state.acc = zeros(3,1);
    elseif(t == 0)
        desired_state.pos = waypoints0(:,1);
        desired_state.vel = zeros(3,1);
        desired_state.acc = zeros(3,1);
    else
        t_index = find(traj_time >= t,1) - 1;
        t = t - traj_time(t_index);
        scale = t/d0(t_index);
        desired_state.pos(1,1) = polyval(alpha_x(t_index,:),scale);
        desired_state.pos(2,1) = polyval(alpha_y(t_index,:),scale);
        desired_state.pos(3,1) = polyval(alpha_z(t_index,:),scale);
                
        desired_state.vel(1,1) = polyval(polyder(alpha_x(t_index,:)),scale)/d0(t_index);
        desired_state.vel(2,1) = polyval(polyder(alpha_y(t_index,:)),scale)/d0(t_index);
        desired_state.vel(3,1) = polyval(polyder(alpha_z(t_index,:)),scale)/d0(t_index);
        
        desired_state.acc(1,1) = polyval(polyder(polyder(alpha_x(t_index,:))),scale)/d0(t_index)^2;
        desired_state.acc(2,1) = polyval(polyder(polyder(alpha_y(t_index,:))),scale)/d0(t_index)^2;
        desired_state.acc(3,1) = polyval(polyder(polyder(alpha_z(t_index,:))),scale)/d0(t_index)^2;
    end
    desired_state.yaw = 0;
    desired_state.yawdot = 0;
end
end

