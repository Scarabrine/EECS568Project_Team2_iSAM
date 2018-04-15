function [] = iSAMRelinearization_pre()
% iSAM relinearization for preprocessed data

global State
global Param

% 1. solve A * delta =b
delta = State.iSAM.R \ State.iSAM.b;
% 2. new_linearization_point = delta + prev_linearization_point
State.iSAM.rlin = State.iSAM.rlin + delta(State.iSAM.rM)';
State.iSAM.llin = State.iSAM.llin + delta(State.iSAM.lM)';
for ind = 1:1:State.iSAM.nR
    State.iSAM.rlin(ind*3) = minimizedAngle(State.iSAM.rlin(ind*3));
end
% 3. relinearize all poses and lardmarks
State.iSAM.R = sparse(size(State.iSAM.R,1), size(State.iSAM.R,2));
State.iSAM.b = sparse(size(State.iSAM.b,1), size(State.iSAM.b,2));
State.iSAM.rM = 1:State.iSAM.nR*3;
State.iSAM.lM = State.iSAM.nR*3 + (1:State.iSAM.nL*2);
pre_row_id = 0;
for ind = 1:1:State.iSAM.nR
    if ind == 1
        Sigma_u = eye(3) * 0.01; %% Anchor sigma, Need to check what's in their implementation
        Lu = pinv(chol(Sigma_u,'lower'));
        G = -eye(3);
        State.iSAM.R((pre_row_id+1):(pre_row_id+3),State.iSAM.rM((ind*3-2):(ind*3))) = Lu * G;
        State.iSAM.b((pre_row_id+1):(pre_row_id+3),1) = Lu * zeros(3,1);
        pre_row_id = pre_row_id +3;
    else
        x1 = State.iSAM.rlin((ind*3-5):(ind*3-3));
        x2 = State.iSAM.rlin((ind*3-2):(ind*3));
        u = State.iSAM.uList(1:3, ind-1);
        c = cos(x1(3));
        s = sin(x1(3));
        dx = x1(1)-x2(1); dy = x1(2) - x2(2);
        G = [c, s, 0;
            -s, c, 0;
             0, 0, 1];
        F = [-c, -s, s*dx-c*dy; 
              s, -c, c*dx+s*dy;
              0,  0, -1];
        % Upate R matrix
        State.iSAM.R((pre_row_id+1):(pre_row_id+3), State.iSAM.rM((ind*3-2):(ind*3))) =  Param.Inf_sqrt_u * G;
        State.iSAM.R((pre_row_id+1):(pre_row_id+3), State.iSAM.rM((ind*3-5):(ind*3-3))) = Param.Inf_sqrt_u * F;
        % Update b vector
        u_hat = [(-dx*c-dy*s);(-dy*c+dx*s);(x2(3)-x1(3))];
        b_cur = u_hat -u;
        b_cur(3) = minimizedAngle(b_cur(3));
        State.iSAM.b((pre_row_id+1):(pre_row_id+3),1) = -Param.Inf_sqrt_u * b_cur;
        pre_row_id = pre_row_id + 3;
    end
end
% relinearize landmark pose
for ind =1:1:State.iSAM.nL
    obsList = State.iSAM.iL{ind};
    for ind_obs = 1:1:size(obsList,2)
        landmark_index = State.iSAM.lM((ind*2-1):(ind*2));
        obs_index = obsList(ind_obs);
        robot_pose_id = State.iSAM.zList(3,obs_index);
        robot_index = State.iSAM.rM((robot_pose_id*3-2):(robot_pose_id*3));
        cur_robot_lin = State.iSAM.rlin((robot_pose_id*3-2):(robot_pose_id*3));
        cur_land_lin = State.iSAM.llin((ind*2-1):(ind*2));
        
        % compute H, J
        c = cos(cur_robot_lin(3)); s = sin(cur_robot_lin(3)); 
        dx = cur_land_lin(1) - cur_robot_lin(1); dy = cur_land_lin(2) - cur_robot_lin(2);
        H = [-c, -s, -s*dx + c*dy;
              s, -c, -c*dx - s*dy];
        J = [c, s;
            -s, c];
        % Update R, b
        State.iSAM.R((pre_row_id+1):(pre_row_id+2), landmark_index) = Param.Inf_sqrt_z * J;
        State.iSAM.R((pre_row_id+1):(pre_row_id+2), robot_index) = Param.Inf_sqrt_z * H;
        z = State.iSAM.zList(1:2,obs_index);
        z_hat = [c, s; -s, c] * [dx;dy];
        b_cur = z_hat - z;
        State.iSAM.b((pre_row_id+1):(pre_row_id+2),1) = -Param.Inf_sqrt_z * b_cur;
        pre_row_id = pre_row_id + 2;
    end
end

% 4. reorder
p = colamd(State.iSAM.R);

for ind = 1:1:length(State.iSAM.rM)
    State.iSAM.rM(ind) = ( find(p == State.iSAM.rM(ind)));
end

for ind = 1:1:length(State.iSAM.lM)
    State.iSAM.lM(ind) = ( find(p == State.iSAM.lM(ind)));
end
% 5. QR decomposition
[Q, R] = qr(State.iSAM.R(:,p));
State.iSAM.R = sparse(R);
State.iSAM.b = sparse(Q') * sparse(State.iSAM.b);
% [m,n] = size(State.iSAM.R);
% givenRotationUpdateRb(1,n,1, m);

end

