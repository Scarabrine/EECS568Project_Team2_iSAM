function [] = iSAMInit()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global Param;
global State;

State.iSAM.rlin = Param.initialStateMean';
State.iSAM.rM = 1:3;
Sigma_u = eye(3)*0.01; %% TODO: set this in the Params
Lu = chol(Sigma_u,'lower');
G = -eye(3);
State.iSAM.R = sparse(Lu \ G);
State.iSAM.Lu = [Lu];
State.iSAM.b = sparse(Lu \ zeros(3,1));
State.iSAM.nR = 1;

end

