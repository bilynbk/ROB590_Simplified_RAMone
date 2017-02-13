function tau = groundController(x,phase,t_prev_stance,k_des,dx_des)
n = size(x,1);
tau = zeros(n/2,1);

%% parameters
param = yumingParameters();
sysParam = param.sysParam;

%% %%%%%%%%%%%%%%% Stance foot %%%%%%%%%%%%%%%
%% Virtual force
F = zeros(3,1);

if (phase == 2) || (phase == 3)
    theta = ThetaR(x(1:n/2),sysParam);
    L_sp = SpringLengthR(x(1:n/2),sysParam);
elseif (phase == 5) || (phase == 6)
    theta = ThetaL(x(1:n/2),sysParam);
    L_sp = SpringLengthL(x(1:n/2),sysParam);
end

k = param.k;
if (phase == 3) || (phase == 6)
    k = k_des;
end
F(1) = -k*(param.L_sp0-L_sp)*sin(theta);
F(2) =  k*(param.L_sp0-L_sp)*cos(theta);

maxForce = 3000;%300;
if (F(1)^2+F(2)^2)>maxForce^2
    F(1) = F(1)*maxForce/(F(1)^2+F(2)^2)^0.5;
    F(2) = F(2)*maxForce/(F(1)^2+F(2)^2)^0.5;
end

%% Body balance (stance foot hip joint)
if phase == 2 || phase == 5
    tar_angle = 0; 
    % PD controller parameters
    kp = 100;    % 10   % 200
    kd = 4;      % 0.5  % 5
    max_f = 1000;    % maximum torque that can be applied
    % PD controller for desired phi.
    err = x(3) - tar_angle;
    derr =  x(3+n/2);     
    tau_balance = -1*(kp*err + kd*derr);
elseif phase == 3 || phase == 6
    % I found the robot lean forward too much with tar_vel = 0.
    tar_vel = 1.2*dx_des + 0.3;
                % can change dx_des to x(6), so there is no delay.
                % if dx_des = 1, then tar_vel should be 1.5
                % if dx_des = 0, then tar_vel should be 0.3
    % Failed try:
    % tar_vel = -((x(8)+x(9))*param.J2+(x(8)+x(9)+x(10))*param.J3)/param.J1;
    
    % P controller parameters
    kp = 10;%*50;       % 2
    % kd = 0.2;   % 0.2
    max_f = 1000;    % maximum torque that can be applied
    % P controller for desired angular velocity.
    err = x(3+n/2) - tar_vel;
    % derr =  x(5)-x(6);
    tau_balance = -kp*err;% - kd*derr;
end
if tau_balance > max_f
    tau_balance = max_f;
elseif tau_balance < -max_f
    tau_balance = -max_f;
end
% Assignment
F(3) = tau_balance;

%% Virtual Force Conversion (when the virtual spring is between the body and the foot)
if phase == 2 || phase == 3
    phi_f = -(x(3)+x(4)+x(5));
    JT = J_virtual_force(phi_f,x(5),x(4),param.lH,param.lL2,param.lL3);
    tau_kh = JT*[F(2);F(3)];
    % Hip joint
    tau(4) = tau_kh(2);
    % Knee joint
    tau(5) = tau_kh(1);
elseif phase == 5 || phase == 6
    phi_f = -(x(3)+x(6)+x(7));
    JT = J_virtual_force(phi_f,x(7),x(6),param.lH,param.lL2,param.lL3);
    tau_kh = JT*[F(2);F(3)];
    % Hip joint
    tau(6) = tau_kh(2);
    % Knee joint
    tau(7) = tau_kh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% %%%%%%%%%%%%%%% Swing foot %%%%%%%%%%%%%%%
%% Next flight phase angle (Get the desired touch-down angle)
% Position controller parameters
kp_pos = param.k_f(1);
kd_pos = param.k_f(2);
% Position controller (PD control)
dx_des = -kp_pos*(x(1)-param.target_pos) - kd_pos*x(8); % target_pos is fixed.
if dx_des>param.max_dx_des
    dx_des = param.max_dx_des;
elseif dx_des<-param.max_dx_des
    dx_des = -param.max_dx_des;
end

% dx_des = 1; % for tuning kp_rai
% dx_des = 0.1; % for tuning kp_rai
% dx_des = 0; % for testing

% Raibert style controller parameters
kp_rai = param.k_f(3);
max_theta_tar = 50*pi/180;
% Raibert style controller
x_des = x(1+n/2)*t_prev_stance/2 + kp_rai*(x(1+n/2)-dx_des);
theta_tar = asin(x_des/param.L_sp0);
if theta_tar > max_theta_tar
    theta_tar = max_theta_tar;
elseif theta_tar < -max_theta_tar
    theta_tar = -max_theta_tar;
end

% theta_tar = 1; % for tuning the next PD controller
% theta_tar = 0; % for debugging

%% Swing Foot: Hip joint (Use "theta and desired theta" instead of alpha and desired alpha)
% Derive current theta (spring is between CoG and foot!)
if phase == 2 || phase == 3 
    theta = ThetaL(x(1:n/2),sysParam);
    d_theta = dThetaL(x,sysParam);
elseif phase == 5 || phase == 6
    theta = ThetaR(x(1:n/2),sysParam);
    d_theta = dThetaR(x,sysParam);
end
% PD controller parameters
kp = 120;    % 120 % 300
kd = 4.8;    % 4.8 % 8
max_f = 1000;    % maximum torque that can be applied
% PD controller for desired phi.
err = theta - theta_tar;
derr =  d_theta;     
        %%% TODO: theta_target is dynamic, so I should modify derr.
tau_hip = -kp*err - kd*derr;
if tau_hip > max_f
    tau_hip = max_f;
elseif tau_hip < -max_f
    tau_hip = -max_f;
end

% Assignment
if phase == 2 || phase == 3 
    tau(6) = tau_hip;
elseif phase == 5 || phase == 6 
    tau(4) = tau_hip;
end

%% Swing Foot: Knee joint
% PD controller parameters
kp = 100;    
kd = 0.9;   
max_f = 1000;    % maximum torque that can be applied
% PD controller for desired phi.
if phase == 2 || phase == 3
    err = x(7) - param.beta_eq*5;
    derr =  x(7+n/2);   
elseif phase == 5 || phase == 6
    err = x(5) - param.beta_eq*5;
    derr =  x(5+n/2);      
end
tau_knee = -kp*err - kd*derr;
if tau_knee > max_f
    tau_knee = max_f;
elseif tau_knee < -max_f
    tau_knee = -max_f;
end

% Assignment
if phase == 2 || phase == 3
    tau(7) = tau_knee;
elseif phase == 5 || phase == 6
    tau(5) = tau_knee;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Limit
tau_max = param.tau_max;
if param.torque_limit_flag
    for i = 4:7
        if abs(tau(i))>tau_max
            tau(i) = sign(tau(i))*tau_max;
        end
    end
end

%% Testing

end