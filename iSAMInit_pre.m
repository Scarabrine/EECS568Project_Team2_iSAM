function [] = iSAMInit_pre()
% Initialize iSAM using preprocess data

global Param;
global State;

State.iSAM.rlin = Param.initialStateMean';
State.iSAM.rM = 1:3;
Sigma_u = eye(3) * 0.01; %% Anchor sigma, Need to check what's in their implementation
Lu = pinv(chol(Sigma_u,'lower'));
G = -eye(3);
State.iSAM.R = Lu * G;
State.iSAM.b = Lu * zeros(3,1);
State.iSAM.nR = 1;
State.iSAM.max_node = 0; %% the maximum index of the nodes.
State.iSAM.rNode = 0;
end

