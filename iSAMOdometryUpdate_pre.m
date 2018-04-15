function [] = iSAMOdometryUpdate_pre(u,nodes)
% iSAMOdometryUpdate for preprocessed vp data
    global Param;
    global State;
    %% Check if the nodes are correct, not useful in the future, commented to speed up
    
%     if nodes(2) == (State.iSAM.max_node+1) && nodes(1) == State.iSAM.rNode(end)
%     else
%         disp('FUCK! Odometry nodes do not match!!!!')
%     end
    
    x1 = State.iSAM.rlin(State.iSAM.nR*3-2:State.iSAM.nR*3);
    c = cos(x1(3));
    s = sin(x1(3));
    x2 = [u(1)*c-u(2)*s+x1(1), u(1)*s+u(2)*c+x1(2), minimizedAngle(u(3)+x1(3))];
    dx = x1(1)-x2(1); dy = x1(2) - x2(2);
    G = [c, s, 0;
        -s, c, 0;
         0, 0, 1];
    F = [-c, -s, s*dx-c*dy; 
          s, -c, c*dx+s*dy;
          0,  0, -1];
      
    % Update R matrix
    State.iSAM.R = blkdiag(State.iSAM.R, Param.Inf_sqrt_u * G);
    State.iSAM.R(end-2:end, State.iSAM.rM(end-2:end)) = Param.Inf_sqrt_u * F;
    % Update rM, rlin, rNode, nR, max_node
    state_dim = size(State.iSAM.R,2);
    State.iSAM.rM = [State.iSAM.rM, state_dim-2:state_dim];
    State.iSAM.rlin = [State.iSAM.rlin, x2];
    State.iSAM.rNode = [State.iSAM.rNode, nodes(2)];
    State.iSAM.nR = State.iSAM.nR + 1;
    State.iSAM.max_node = State.iSAM.max_node + 1;
    % Update b vector
    u_hat = [(-dx*c-dy*s);(-dy*c+dx*s);(x2(3)-x1(3))];
    b_cur = u_hat -u';
    b_cur(3) = minimizedAngle(b_cur(3));
    State.iSAM.b = [State.iSAM.b; -Param.Inf_sqrt_u * b_cur];
    % Keep R to be up triangle, given rotation is not quick enough!!!!!
    [Q, R] = qr(State.iSAM.R);
    State.iSAM.R = sparse(R);
    State.iSAM.b = sparse(Q') * sparse(State.iSAM.b);
    %% These code is for given rotation, not fast enough, not check in this file.
%     [m,n] = size(State.iSAM.R);
%     k_start = State.iSAM.rM(end-5);% start of F colomn dimension
%     k_end = n;
%     i_start = m - 2;
%     i_end = m;
%     givenRotationUpdateRb(k_start, k_end, i_start, i_end);
    State.iSAM.uList = [State.iSAM.uList,u'];
  
end

