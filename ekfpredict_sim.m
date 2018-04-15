function ekfpredict_sim(u)
% EKF-SLAM prediction for simulator process model

global Param;
global State;

% Get mu, sigma for the augmented state vector for the previous step
mu_t_1 = State.Ekf.mu; % Get mu_{t-1}
Sigma_t_1 = State.Ekf.Sigma; % Get Sigma_{t-1}

% Get number of landmarks and state dimension
num_landmark = State.Ekf.nL; % Need to be care reset this in the correction step
num_state = length(mu_t_1); % state vector dimension

% Set Fx matrix
Fx = [eye(3) zeros(3, num_landmark * 2)];
% Get control input and robot pose
drot1 = u(1); dtrans = u(2); drot2 = u(3);
x = mu_t_1(1); y = mu_t_1(2); theta = mu_t_1(3);
% Process vector for the robot
g = [dtrans * cos(theta + drot1); dtrans * sin(theta + drot2); drot1 + drot2];

% Compute the mu_t
mu_t = mu_t_1 + Fx' * g;
mu_t(3) = minimizedAngle(mu_t(3)); % normalize angle after prediction step
% compute Gt
Gtx = [0 0 -dtrans * sin(theta + drot1);
       0 0 dtrans * cos(theta + drot2);
       0 0 0]; % This Gtx is G_t^x - I in the writeup.
Gt = eye(num_state) + Fx' * Gtx * Fx;

% compute Vt
Vt = [-dtrans * sin(theta + drot1) cos(theta + drot1) 0;
       dtrans * cos(theta + drot1) sin(theta + drot1) 0;
       1 0 1];
% compute Mt
a1 = Param.alphas(1); % parameter alpha1
a2 = Param.alphas(2); % parameter alpha2
a3 = Param.alphas(3); % parameter alpha3
a4 = Param.alphas(4); % parameter alpha4
Mt = [a1*drot1^2 + a2 * dtrans^2 0                                    0;
      0                          a3*dtrans^2 + a4*(drot1^2 + drot2^2) 0;
      0                          0                                    a1*drot2^2 + a2*dtrans^2];
% compute Rtx
Rtx = Vt * Mt * Vt';

% compute Sigma_t
Sigma_t = Gt * Sigma_t_1 * Gt' + Fx' * Rtx * Fx;

% Set mu_t and sigma_t back to the State global variable.
State.Ekf.mu = mu_t;
State.Ekf.Sigma = Sigma_t;
end
