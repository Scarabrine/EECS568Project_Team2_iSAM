function ekf_pre_update (z, nodes)

global Param;
global State;

switch lower(Param.dataAssociation)
    case 'known'
        Li = da_known_pre(z, nodes);
    otherwise
        error('not test for this data association method!!!!')
end

nPoseState = 3;
state_size = length(State.Ekf.mu);

if(Li>0)

    F = zeros(5,state_size);
    F(1:3,1:3) = eye(3);
    F(4:5, 2*Li+2:2*Li+3) = eye(2);
    
    % expected observation
    c = cos(State.Ekf.mu(3)); s = sin(State.Ekf.mu(3)); 
    dx = State.Ekf.mu(2*Li+2) - State.Ekf.mu(1); 
    dy = State.Ekf.mu(2*Li+3) - State.Ekf.mu(2);

    z_hat = [c, s; -s, c] * [dx;dy];    
    
    H = Jacobian_observation(Li) * F;
    
    % Kalman Gain
    C = (H*State.Ekf.Sigma*H'+Param.R);
    C = C\eye(size(C));
    K = State.Ekf.Sigma * H' * C;
    % calculate change in mu, sigma
    obs_error = z(1:2)'-z_hat(1:2);
    obs_error(2) = minimizedAngle(obs_error(2));
    
    State.Ekf.mu = State.Ekf.mu + K * obs_error;
    State.Ekf.mu(3) = minimizedAngle(State.Ekf.mu(3));
    State.Ekf.Sigma = (eye(size(State.Ekf.Sigma)) - K*H ) * State.Ekf.Sigma;
end
if(Li==0)
    initialize_new_landmark(z, nodes);
end
% plotting estimate ellipses
% plotcov2d(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma(1:2,1:2), 'k', 1, 'k', 0.3, 3);
% for i=4:2:length(State.Ekf.mu)
%     plotcov2d(State.Ekf.mu(i), State.Ekf.mu(i+1), State.Ekf.Sigma(i:i+1,i:i+1), 'r', 1, 'r', 0.1, 3);
% end