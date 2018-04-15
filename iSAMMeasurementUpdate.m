function [] = iSAMMeasurementUpdate(z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global Param;
global State;

switch lower(Param.dataAssociation)
    case 'known'
        [Li,measurements] = da_known(z, State.iSAM.sL);
    case 'nndg'
        [Li,measurements]=da_nn_dg_sim_isam(z);
    case 'nn'
        [Li,measurements] = da_nn_isam(z, Param.R);
    case 'jcbb'
        [Li,measurements] = da_jcbb(z, Param.R);
    otherwise
        error('unrecognized data association method: "%s"', Param.dataAssociation);
end

%% Compute Lz
Sigma_z = Param.R;
Lz = chol(Sigma_z, 'lower');


for ind =1:1:size(z,2)
    % check if this observation/landmark is new or not
    if Li(ind) == 0
        %% This landmark is new, add new landmark
        Li(ind) = State.iSAM.nL+1;
        State.iSAM.nL = State.iSAM.nL+1;
        State.iSAM.sL = [State.iSAM.sL, z(3,ind)];%% Should be the most frequent signature.
        %% loop all the measurements observing this new landmark, compute the landmark linearization point
        cur_land_lin = [0,0];
        for ind_obs = 1:1:size(measurements{ind},2)
            robot_id = measurements{ind}(1,ind_obs);
            obs = measurements{ind}(2:4, ind_obs);
            cur_rob_lin = State.iSAM.rlin((3*robot_id-2):3*robot_id);
%             cur_rob_lin = measurements{ind}(5:7,ind_obs);

            r = obs(1);
            phi = obs(2);
            cur_land_lin = cur_land_lin + [cur_rob_lin(1) + r*cos(phi + cur_rob_lin(3)), cur_rob_lin(2) + r*sin(phi + cur_rob_lin(3))];
        end
        %% new landmark linearization point is the mean of those cur_land_lin.
        cur_land_lin = cur_land_lin / size(measurements{ind},2);
        State.iSAM.llin = [State.iSAM.llin, cur_land_lin];
        State.iSAM.lM = [State.iSAM.lM, size(State.iSAM.R,2)+1:size(State.iSAM.R,2)+2];
        State.iSAM.iL{State.iSAM.nL} = [];
        for ind_obs = 1:1:size(measurements{ind},2)
            % get current robot and landmark linearization point
            robot_id = measurements{ind}(1,ind_obs);
            cur_rob_lin = State.iSAM.rlin((3*robot_id-2):3*robot_id);
%             cur_rob_lin = measurements{ind}(5:7,ind_obs);

            cur_obs = measurements{ind}(2:4, ind_obs);

            % compute H,J
            delta = cur_land_lin - cur_rob_lin(1:2);
            q = delta * delta';
            H = 1 / q * [-sqrt(q)*delta(1), -sqrt(q)*delta(2), 0;
                                 delta(2),  -delta(1), -q ]; %% robot pose
            J = 1 / q * [  sqrt(q)*delta(1), sqrt(q)*delta(2);
                           -delta(2),        delta(1)        ]; %% landmark pose
            % update R,b
%             State.iSAM.R = blkdiag(State.iSAM.R, Lz \ J);
            State.iSAM.R(end+1:end+2, State.iSAM.lM(end-1:end)) = Lz \ J;
            State.iSAM.R(end-1:end, State.iSAM.rM((robot_id*3-2):robot_id*3)) = Lz \ H;

            z_hat = [sqrt(q); atan2(delta(2), delta(1)) - cur_rob_lin(3)];
            z_hat(2) = minimizedAngle(z_hat(2));
            b_cur = cur_obs(1:2) - z_hat;
            b_cur(2) = minimizedAngle(b_cur(2)); 
            State.iSAM.b = [State.iSAM.b; Lz \ b_cur];
            % given rotation, update R,b
            [m,n] = size(State.iSAM.R);
            k_start = State.iSAM.rM(robot_id*3-2);
            k_end = n;
            i_start = m-1;
            i_end = m;
    %         givenRotationUpdateRb(k_start, k_end, i_start, i_end);
            [Q, R] = qr(State.iSAM.R);
            State.iSAM.R = sparse(R);
            State.iSAM.b = sparse(Q') * sparse(State.iSAM.b);

%             State.iSAM.Lz = blkdiag(State.iSAM.Lz, Lz);
            State.iSAM.zList = [State.iSAM.zList, [obs(1:2);robot_id]];
            State.iSAM.iL{State.iSAM.nL} = [State.iSAM.iL{State.iSAM.nL}, size(State.iSAM.zList,2)];
        end
    elseif Li(ind) >0
        %% This landmark has been observed.
        landmark_id = Li(ind);
        landmark_index = State.iSAM.lM((landmark_id*2-1):landmark_id*2);
        robot_index = State.iSAM.rM(end-2:end);
        
        cur_rob_lin = State.iSAM.rlin(end-2:end);
        cur_land_lin = State.iSAM.llin((landmark_id*2-1):landmark_id*2);
        
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
        
        z_hat = [sqrt(q); atan2(delta(2), delta(1)) - cur_rob_lin(3)];
        z_hat(2) = minimizedAngle(z_hat(2));
        b_cur = z(1:2,ind) - z_hat;
        b_cur(2) = minimizedAngle(b_cur(2)); 
        State.iSAM.b = [State.iSAM.b; Lz \ b_cur];
        % given rotation, update R,b
        [m,n] = size(State.iSAM.R);
        k_start = min([robot_index, landmark_index]);
        k_end = n;
        i_start = m - 1;
        i_end = m;
%         givenRotationUpdateRb(k_start, k_end, i_start, i_end);
        [Q, R] = qr(State.iSAM.R);
        State.iSAM.R = sparse(R);
        State.iSAM.b = sparse(Q') * sparse(State.iSAM.b);
        
        State.iSAM.zList = [State.iSAM.zList, [z(1:2,ind);State.iSAM.nR]];
        State.iSAM.iL{landmark_id} = [State.iSAM.iL{landmark_id}, size(State.iSAM.zList,2)];
    else
        % Li(ind) == -1, do nothing.
    end
    
    
end



end

