function ekfupdate(z)
% EKF-SLAM update step for both simulator and Victoria Park data set

global Param;
global State;


% returns state vector indices pairing observations with landmarks
switch lower(Param.dataAssociation)
    case 'known'
        Li = da_known(z,State.Ekf.sL); 
    case 'nn'
        Li = da_nn(z(1:2,:), Param.R);
    case 'nndg'
        [Li,measurements] = da_nn_dg_sim_ekf(z);
        for i=1:size(z,2)
            if Li(i)==0
                z(:,i)=measurements{i};
            end
        end
    case 'jcbb'
        [Li,~] = da_jcbb(z(1:2,:), Param.R);
    otherwise
        error('unrecognized data association method: "%s"', Param.dataAssociation);
end


switch lower(Param.updateMethod)
    case 'batch'
        % get Qt, sensing noise covariance matrix
        Q_t_i = Param.R;
        H_t = [];
        delta_t = [];
        Q_t = [];
        %% Check if there are new landmarks.
         if Li(ind)~=-1
        for ind = 1:1:size(z,2)
            %% Check if the landmark is new or not
           
            if Li(ind) == 0
               %% This landmark is new, add new landmark
               initialize_new_landmark_ekf(z(:,ind), Param.R);
               Li(ind) = State.Ekf.nL;
            end
            State.Ekf.iL{Li(ind)} = [State.Ekf.iL{Li(ind)}, State.Ekf.t];
        end
        
        for ind = 1:1:size(z,2)
            % get robot pose
            pos_r = State.Ekf.mu(1:2);
            theta = State.Ekf.mu(3);
            landmark_id = Li(ind);
            pos_landmark = State.Ekf.mu(State.Ekf.iM((landmark_id*2-1):landmark_id*2));
            delta = pos_landmark - pos_r;
            q = delta'*delta;
            % compute the estimated observation
            z_hat_t_i = [sqrt(q);atan2(delta(2),delta(1)) - theta];
            z_hat_t_i(2) = minimizedAngle(z_hat_t_i(2));
            % compute the projection matrix F_x_j
            F_x_j = zeros(5,length(State.Ekf.mu));
            F_x_j(1:3,1:3) = eye(3);
            F_x_j(4:5,State.Ekf.iM((landmark_id*2-1):landmark_id*2)) = eye(2);
            % compute H_t_i
            H_t_i = 1 / q * [-sqrt(q)*delta(1) -sqrt(q)*delta(2) 0  sqrt(q)*delta(1) sqrt(q)*delta(2);
                             delta(2)          -delta(1)         -q -delta(2)        delta(1)        ] * F_x_j;
                         
            % Stack H_t
            H_t = [H_t;H_t_i];
                
            % Stack delta_t;
            delta_z = z(1:2,ind) - z_hat_t_i;
            delta_z(2) = minimizedAngle(delta_z(2));
            delta_t = [delta_t;delta_z];
            
            % Stack Q_t
            Q_t = blkdiag(Q_t,Q_t_i);
        end
        % Batch update.
        try
            S_t = H_t * State.Ekf.Sigma * H_t' + Q_t;
        catch
            disp('problem')
        end
        K_t = State.Ekf.Sigma * H_t'*pinv(S_t);
        State.Ekf.mu = State.Ekf.mu + K_t*delta_t;
        State.Ekf.Sigma = (eye(size(State.Ekf.Sigma,1)) - K_t*H_t)*State.Ekf.Sigma;
        end
    case 'seq' % TODO: need to check for seq, the initialize landmark should be in the loop of check observations. not in the da.
        % get Qt, sensing noise covariance matrix
        Q_t_i = Param.R;
        for ind = 1:1:size(z,2)
            %% Check if the landmark is new or not
            if Li(ind)~=-1
            if Li(ind) == 0
               %% This landmark is new, add new landmark
               initialize_new_landmark_ekf(z(:,ind), Param.R);
               Li(ind) = State.Ekf.nL;
            end

            % get robot pose
            pos_r = State.Ekf.mu(1:2);
            theta = State.Ekf.mu(3);
            
            landmark_id = Li(ind);
            pos_landmark = State.Ekf.mu(State.Ekf.iM((landmark_id*2-1):landmark_id*2));
            delta = pos_landmark - pos_r;
            q = delta'*delta;
            % compute the estimated observation
            z_hat_t_i = [sqrt(q);atan2(delta(2),delta(1)) - theta];
            z_hat_t_i(2) = minimizedAngle(z_hat_t_i(2));
            % compute the projection matrix F_x_j
            F_x_j = zeros(5,length(State.Ekf.mu));
            F_x_j(1:3,1:3) = eye(3);
            F_x_j(4:5,State.Ekf.iM((landmark_id*2-1):landmark_id*2)) = eye(2);
            % compute H_t_i
            H_t_i = 1 / q * [-sqrt(q)*delta(1) -sqrt(q)*delta(2) 0  sqrt(q)*delta(1) sqrt(q)*delta(2);
                             delta(2)          -delta(1)         -q -delta(2)        delta(1)        ] * F_x_j;
            % compute the karman gain K_t_i
            K_t_i = State.Ekf.Sigma * H_t_i' * pinv(H_t_i * State.Ekf.Sigma * H_t_i' + Q_t_i);
            % update mu and sigma
            delta_z = z(1:2,ind) - z_hat_t_i;
            delta_z(2) = minimizedAngle(delta_z(2));
            State.Ekf.mu = State.Ekf.mu + K_t_i * delta_z;
            State.Ekf.Sigma = (eye(size(State.Ekf.Sigma,1)) - K_t_i * H_t_i) * State.Ekf.Sigma;
            if State.Ekf.t == 1 && ind == 1
               %% Get the first landmark sigma
               State.Ekf.mu_first = State.Ekf.mu;
               State.Ekf.Sigma_first = State.Ekf.Sigma;
            end
            end
        end
    otherwise
        error('unrecognized update method: "%s"', Param.updateMethod);
end









end