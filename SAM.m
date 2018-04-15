function [] = SAM()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global State
global Param
% 3. relinearize all poses and landmarks
State.iSAM.R = [];
State.iSAM.b = [];
State.iSAM.Lu = [];
State.iSAM.Lz = [];
State.iSAM.rM = [];
State.iSAM.lM = [];
% relinearize robot poses
for ind =1:1:State.iSAM.nR
    if ind == 1
        G = -eye(3);

        Sigma_u = eye(3)*0.01;
        Lu = chol(Sigma_u,'lower');
        State.iSAM.R = Lu \ G;
        State.iSAM.b = Lu \ zeros(3,1); %% Need double check this problem.
        State.iSAM.Lu = Lu;
        State.iSAM.rM = 1:3;
    else
        % get current linearization point and control
        u = State.iSAM.uList(:,ind-1);
        drot1 = u(1); dtrans = u(2); drot2 = u(3);
        cur_rob_lin = State.iSAM.rlin((ind*3-2):ind*3);
        prev_rob_lin = State.iSAM.rlin((ind*3-5):(ind*3-3));
        x_prev = prev_rob_lin(1);
        y_prev = prev_rob_lin(2);
        theta_prev = prev_rob_lin(3);

        % 1. compute Lu
        a1 = Param.alphas(1); % parameter alpha1
        a2 = Param.alphas(2); % parameter alpha2
        a3 = Param.alphas(3); % parameter alpha3
        a4 = Param.alphas(4); % parameter alpha4
        M = [ a1*drot1^2 + a2 * dtrans^2 0                                    0;
              0                          a3*dtrans^2 + a4*(drot1^2 + drot2^2) 0;
              0                          0                                    a1*drot2^2 + a2*dtrans^2];
        V = [-dtrans * sin(theta_prev + drot1) cos(theta_prev + drot1) 0;
           dtrans * cos(theta_prev + drot1) sin(theta_prev + drot1) 0;
           1 0 1];
%         Sigma_u = V*M*V' + eye(3) * 0.00001;
        Sigma_u = eye(3)*0.01;
        Lu = chol(Sigma_u,'lower');
        State.iSAM.Lu = blkdiag(State.iSAM.Lu, Lu);
        % compute F,G
        G = -eye(3);
        F = [1 0 -dtrans * sin(theta_prev + drot1);
               0 1 dtrans * cos(theta_prev + drot2);
               0 0 1];
        % update A,b
        State.iSAM.R = blkdiag(State.iSAM.R, Lu \ G);
        State.iSAM.R(end-2:end, State.iSAM.rM(end-2:end)) = Lu \ F;
        state_dim = size(State.iSAM.R,2);
        State.iSAM.rM = [State.iSAM.rM, state_dim-2:state_dim];
        cur_b = cur_rob_lin' - prediction(prev_rob_lin',u);
        cur_b(3) = minimizedAngle(cur_b(3));
        State.iSAM.b = [State.iSAM.b; Lu \ cur_b];
        
    end
end

% relinearize landmark poses
Sigma_z = Param.R;
Lz = chol(Sigma_z, 'lower');
for ind =1:1:State.iSAM.nL
    obsList = State.iSAM.iL{ind};
    State.iSAM.lM = [State.iSAM.lM, size(State.iSAM.R,2)+1:size(State.iSAM.R,2)+2];
    State.iSAM.Lz = blkdiag(State.iSAM.Lz, Lz);
    for ind_obs = 1:1:size(obsList,2)
        landmark_index = State.iSAM.lM(end-1:end);
        obs_index = obsList(ind_obs);
        robot_pose_id = State.iSAM.zList(3,obs_index);
        robot_index = State.iSAM.rM(robot_pose_id*3-2:robot_pose_id*3);
        cur_rob_lin = State.iSAM.rlin(robot_pose_id*3-2:robot_pose_id*3);
        cur_land_lin = State.iSAM.llin(ind*2-1:ind*2);
        
        % compute H,J
        delta = cur_land_lin - cur_rob_lin(1:2);
        q = delta * delta';
        H = 1 / q * [-sqrt(q)*delta(1), -sqrt(q)*delta(2), 0;
                             delta(2),  -delta(1), -q ]; %% robot pose
        J = 1 / q * [  sqrt(q)*delta(1), sqrt(q)*delta(2);
                       -delta(2),        delta(1)        ]; %% landmark pose
        
        % update R,b
        State.iSAM.R(end+1:end+2, landmark_index) = Lz \ J;
        State.iSAM.R(end-1:end, robot_index) = Lz \ H;  
        
        z = State.iSAM.zList(1:2,obs_index);
        
        z_hat = [sqrt(q); atan2(delta(2), delta(1)) - cur_rob_lin(3)];
        z_hat(2) = minimizedAngle(z_hat(2));
        b_cur = z - z_hat;
        b_cur(2) = minimizedAngle(b_cur(2)); 
        State.iSAM.b = [State.iSAM.b; Lz \ b_cur];
    end
end

% 4. reorder
% p = colamd(State.iSAM.R);
% 
% for ind = 1:1:length(State.iSAM.rM)
%     State.iSAM.rM(ind) = ( find(p == State.iSAM.rM(ind)));
% end
% 
% for ind = 1:1:length(State.iSAM.lM)
%     State.iSAM.lM(ind) = ( find(p == State.iSAM.lM(ind)));
% end

% 5. QR decomposition
% [Q, State.iSAM.R] = qr(State.iSAM.R(:,p));
% State.iSAM.b = Q' * State.iSAM.b;
% [m,n] = size(State.iSAM.R);
% givenRotationUpdateRb(1,n,1, m);


% 1. solve A * delta = b
delta = State.iSAM.R \ State.iSAM.b;
State.iSAM.rlin = State.iSAM.rlin + delta(State.iSAM.rM)';
State.iSAM.llin = State.iSAM.llin + delta(State.iSAM.lM)';
for ind = 1:1:State.iSAM.nR
    State.iSAM.rlin(ind*3) = minimizedAngle(State.iSAM.rlin(ind*3));
end

end

