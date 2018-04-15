function [] = iSAMOdometryUpdate_vp(u,dt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global Param;
global State;


% Get model parameters
L = Param.L;
H = Param.H;
a = Param.a;
b = Param.b;
Qu = Param.Qu;
Qf = Param.Qf;


% 0. get control input, previous linearization point
ve=u(1); alpha = u(2);
vc = ve/(1-tan(alpha)*H/L); % Note that as the noise model if given on vc and alpha, thus the Jacobian is all resprect to vc not ve.
x_prev = State.iSAM.rlin(end-2);
y_prev = State.iSAM.rlin(end-1);
theta_prev = State.iSAM.rlin(end);


% 1. compute Lu

Vt = [dt*(cos(theta_prev) - 1/L*tan(alpha) * (a*sin(theta_prev) + b*cos(theta_prev))), -dt*vc/L/(cos(alpha))^2 * (a*sin(theta_prev) + b*cos(theta_prev));
      dt*(sin(theta_prev) + 1/L*tan(alpha) * (a*cos(theta_prev) - b*sin(theta_prev))),  dt*vc/L/(cos(alpha))^2 * (a*cos(theta_prev) - b*sin(theta_prev));
      dt/L*tan(alpha),                                             dt*vc/L/(cos(alpha))^2];

Rtx = Vt * Qu * Vt';

Sigma_u = Rtx + Qf; %%%% Need to test what to use for Sigma_u, might be just Qf, or Qf+Qu for short.

% Sigma_u = eye(3)*0.01;
try
    Lu = chol(Sigma_u,'lower');
catch
    disp('FUCK')
end
% State.iSAM.Lu = blkdiag(State.iSAM.Lu, Lu);
% 1. compute F,G
G = -eye(3);
F = [1 0 dt*(-vc*sin(theta_prev) - vc/L*tan(alpha) * ( a*cos(theta_prev) - b*sin(theta_prev)));
     0 1 dt*( vc*cos(theta_prev) + vc/L*tan(alpha) * (-a*sin(theta_prev) - b*cos(theta_prev)));
     0 0 1];
 
g = [dt * (vc*cos(theta_prev) - vc/L*tan(alpha) * (a*sin(theta_prev) + b*cos(theta_prev)));
     dt * (vc*sin(theta_prev) + vc/L*tan(alpha) * (a*cos(theta_prev) - b*sin(theta_prev)));
     dt * vc / L * tan(alpha)];
   
   
% 2. add F, G to the R
State.iSAM.R = blkdiag(State.iSAM.R, Lu \ G);
State.iSAM.R(end-2:end, State.iSAM.rM(end-2:end)) = Lu \ F;
% 3. update rM, r_lin
state_dim = size(State.iSAM.R,2);
State.iSAM.rM = [State.iSAM.rM, state_dim-2:state_dim];
current_lin = ([x_prev;y_prev;theta_prev] + g)';

State.iSAM.rlin = [State.iSAM.rlin, current_lin(1), current_lin(2), current_lin(3)];
% 4. update b
State.iSAM.b = [State.iSAM.b; Lu \ zeros(3,1)]; %% current_lin - prediction(prev_lin, u)
% 5. do given rotation update R and b
[m,n] = size(State.iSAM.R);
k_start = State.iSAM.rM(end-5);% start of F colomn dimension
k_end = n;
i_start = m - 2;
i_end = m;
% givenRotationUpdateRb(k_start, k_end, i_start, i_end)
[Q,R] = qr(State.iSAM.R);
State.iSAM.R = sparse(R);
State.iSAM.b = sparse(Q')*sparse(State.iSAM.b);


State.iSAM.nR = State.iSAM.nR + 1;
State.iSAM.uList = [State.iSAM.uList, [u;dt]];
end

