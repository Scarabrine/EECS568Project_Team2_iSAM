function [] = iSAMOdometryUpdate_sim(u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global Param;
global State;



% 0. get control input, previous linearization point
drot1 = u(1); dtrans = u(2); drot2 = u(3);
x_prev = State.iSAM.rlin(end-2);
y_prev = State.iSAM.rlin(end-1);
theta_prev = State.iSAM.rlin(end);
% 1. compute Lu
% a1 = Param.alphas(1); % parameter alpha1
% a2 = Param.alphas(2); % parameter alpha2
% a3 = Param.alphas(3); % parameter alpha3
% a4 = Param.alphas(4); % parameter alpha4
% M = [ a1*drot1^2 + a2 * dtrans^2 0                                    0;
%       0                          a3*dtrans^2 + a4*(drot1^2 + drot2^2) 0;
%       0                          0                                    a1*drot2^2 + a2*dtrans^2];
% V = [-dtrans * sin(theta_prev + drot1) cos(theta_prev + drot1) 0;
%    dtrans * cos(theta_prev + drot1) sin(theta_prev + drot1) 0;
%    1 0 1];
% Sigma_u = V*M*V' + eye(3) * 0.00001;
Sigma_u = eye(3)*0.01;
try
    Lu = chol(Sigma_u,'lower');
catch
    disp('FUCK')
end
% State.iSAM.Lu = blkdiag(State.iSAM.Lu, Lu);
% 1. compute F,G
G = -eye(3);
F = [1 0 -dtrans * sin(theta_prev + drot1);
       0 1 dtrans * cos(theta_prev + drot2);
       0 0 1];
   
% 2. add F, G to the R
State.iSAM.R = blkdiag(State.iSAM.R, Lu \ G);
State.iSAM.R(end-2:end, State.iSAM.rM(end-2:end)) = Lu \ F;
% 3. update rM, r_lin
state_dim = size(State.iSAM.R,2);
State.iSAM.rM = [State.iSAM.rM, state_dim-2:state_dim];
current_lin = prediction([x_prev;y_prev;theta_prev],u);
State.iSAM.rlin = [State.iSAM.rlin, current_lin(1), current_lin(2), current_lin(3)];
% 4. update b
State.iSAM.b = [State.iSAM.b; Lu \ zeros(3,1)]; %% current_lin - prediction(prev_lin, u)
% 5. do given rotation update R and b
% [m,n] = size(State.iSAM.R);
% k_start = State.iSAM.rM(end-5);% start of F colomn dimension
% k_end = n;
% i_start = m - 2;
% i_end = m;
% givenRotationUpdateRb(k_start, k_end, i_start, i_end);

[Q, R] = qr(State.iSAM.R);
State.iSAM.R = sparse(R);
State.iSAM.b = sparse(Q') * sparse(State.iSAM.b);

State.iSAM.nR = State.iSAM.nR + 1;
State.iSAM.uList = [State.iSAM.uList, u];
end

