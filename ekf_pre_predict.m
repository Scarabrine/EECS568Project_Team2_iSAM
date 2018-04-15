function ekf_pre_predict(u)

global Param;
global State;
nstatePose = 3;

State.iSAM.max_node = State.iSAM.max_node+1;

% mean update
State.Ekf.mu(1) = State.Ekf.mu(1) + u(1) * cos(State.Ekf.mu(3)) - u(2) * sin(State.Ekf.mu(3));
State.Ekf.mu(2) = State.Ekf.mu(2) + u(1) * sin(State.Ekf.mu(3)) + u(2) * cos(State.Ekf.mu(3));
State.Ekf.mu(3) = State.Ekf.mu(3) + u(3);

State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));

%linearization
[G,V] = Jacobian_motion (u);

% covariance update
State.Ekf.Sigma(1:nstatePose,1:nstatePose) = G * State.Ekf.Sigma(1:nstatePose,1:nstatePose) * G' + V * Param.M * V';