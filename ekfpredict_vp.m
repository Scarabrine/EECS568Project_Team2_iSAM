function ekfpredict_vp(u, dt)
% EKF-SLAM prediction for Victoria Park process model

global Param;
global State;

% Get mu, sigma for the augmented state vector for the previous step
mu_t_1 = State.Ekf.mu; % Get mu_{t-1}
Sigma_t_1 = State.Ekf.Sigma; %Get Sigma_{t-1}

% Get number of landmarks and state dimension
num_landmark = State.Ekf.nL; % Need to be care reset this in the correction step
num_state = length(mu_t_1); % state vector dimension

% Get model parameters
L = Param.L;
H = Param.H;
a = Param.a;
b = Param.b;
Qu = Param.Qu;
Qf = Param.Qf;

% Set Fx matrix
Fx = [eye(3) zeros(3, num_landmark * 2)];
% Get control input and robot pose
ve = u(1); alpha = u(2);
vc = ve/(1-tan(alpha)*H/L); % Note that as the noise model if given on vc and alpha, thus the Jacobian is all resprect to vc not ve.
x = mu_t_1(1); y = mu_t_1(2); phi = mu_t_1(3);
% Process vector for the robot
g = [dt * (vc*cos(phi) - vc/L*tan(alpha) * (a*sin(phi) + b*cos(phi)));
     dt * (vc*sin(phi) + vc/L*tan(alpha) * (a*cos(phi) - b*sin(phi)));
     dt * vc / L * tan(alpha)];
 
% compute the mu_t
mu_t = mu_t_1 + Fx' * g;
mu_t(3) = minimizedAngle(mu_t(3)); % normalize angle after prediction step

% compute Gt
Gtx = [0 0 dt*(-vc*sin(phi) - vc/L*tan(alpha) * ( a*cos(phi) - b*sin(phi)));
       0 0 dt*( vc*cos(phi) + vc/L*tan(alpha) * (-a*sin(phi) - b*cos(phi)));
       0 0 0]; % This Gtx is G_t^x -I in the writeup.
   
Gt = eye(num_state) + Fx' * Gtx * Fx;

% compute Vt
Vt = [dt*(cos(phi) - 1/L*tan(alpha) * (a*sin(phi) + b*cos(phi))), -dt*vc/L/(cos(alpha))^2 * (a*sin(phi) + b*cos(phi));
      dt*(sin(phi) + 1/L*tan(alpha) * (a*cos(phi) - b*sin(phi))),  dt*vc/L/(cos(alpha))^2 * (a*cos(phi) - b*sin(phi));
      dt/L*tan(alpha),                                             dt*vc/L/(cos(alpha))^2];
  
% compute Rtx
Rtx = Vt * Qu * Vt';

% compute Sigma_t
Sigma_t = Gt * Sigma_t_1 * Gt' + Fx' * Rtx * Fx + Fx' * Qf * Fx;

State.Ekf.mu = mu_t;
State.Ekf.Sigma = Sigma_t;

end