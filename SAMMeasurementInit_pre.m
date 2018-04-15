function [] = SAMMeasurementInit_pre(z,nodes)
% SAM Measurement init for preprocessed data.
    global Param;
    global State;
    switch lower(Param.dataAssociation)
        case 'known'
            Li = da_known_pre(z, nodes);
        otherwise
            error('not test for this data association method!!!!')
    end

%% Currently, we only have one observation at once!!!!!!
    if Li == 0
        %% Need to create new landmark
        %% Initialize new landmark location
            robot_index = find(State.iSAM.rNode == nodes(1));
            if length(robot_index) ~= 1
                error('robot node index should be unique in the rNode list');
            end
            cur_robot_lin = State.iSAM.rlin(robot_index*3-2:robot_index*3);
            c = cos(cur_robot_lin(3)); s = sin(cur_robot_lin(3)); x = cur_robot_lin(1); y = cur_robot_lin(2);
            cur_land_lin = [x + c*z(1) - s*z(2), y + s*z(1) + c*z(2)];
            landmark_index = State.iSAM.nL + 1;
        %% Update information
            State.iSAM.llin = [State.iSAM.llin, cur_land_lin];
            State.iSAM.nL = State.iSAM.nL + 1;
            state_dim = size(State.iSAM.R, 2);
            State.iSAM.lM = [State.iSAM.lM, state_dim+1 : state_dim+2];
            State.iSAM.lNode = [State.iSAM.lNode, nodes(2)];
            State.iSAM.max_node = State.iSAM.max_node + 1;
            State.iSAM.iL{landmark_index} = [];
    else
        %% Do not need to create new landmark
            robot_index = find(State.iSAM.rNode == nodes(1));
            if length(robot_index) ~=1
                error('robot node index should be unique in the rNode list');
            end
            cur_robot_lin = State.iSAM.rlin(robot_index*3-2:robot_index*3);
            landmark_index = Li;
            cur_land_lin = State.iSAM.llin(landmark_index*2-1:landmark_index*2);
    end
    % compute H and J
    c = cos(cur_robot_lin(3)); s = sin(cur_robot_lin(3)); 
    dx = cur_land_lin(1) - cur_robot_lin(1); dy = cur_land_lin(2) - cur_robot_lin(2);
    H = [-c, -s, -s*dx + c*dy;
          s, -c, -c*dx - s*dy];
    J = [c, s;
        -s, c];
    %% Update R matrix
    State.iSAM.R(end+1:end+2, State.iSAM.lM(landmark_index*2-1:landmark_index*2)) = Param.Inf_sqrt_z * J;
    State.iSAM.R(end-1:end, State.iSAM.rM(robot_index*3-2:robot_index*3)) = Param.Inf_sqrt_z * H;
    %% Update b matrix
    z_hat = [c, s; -s, c] * [dx;dy];
    b_cur = z_hat - z';
    State.iSAM.b = [State.iSAM.b; -Param.Inf_sqrt_z * b_cur];
    State.iSAM.zList = [State.iSAM.zList, [z';robot_index]];
    State.iSAM.iL{landmark_index} = [State.iSAM.iL{landmark_index}, size(State.iSAM.zList,2)];
    if robot_index ~= State.iSAM.nR
        error('robot_index does not equal to the current robot poses!!!!');
    end
    


end

