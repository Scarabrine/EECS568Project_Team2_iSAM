function [mnew,Sigma]=get_pred_landmark_position(z,current_pose,stateSigma)

global Param
x_k=current_pose;
predSigma=stateSigma;
%get predicted pose
%predicted measurement 
r_k=z(1);
theta_k=z(2);

%estimate new landmark position
%if Param.vp
   %mnew=[r_k*cos(theta_k+x_k(3)-pi/2)+x_k(1);r_k*sin(theta_k+x_k(3)-pi/2)+x_k(2)];
  % g_xk=[1,0,-r_k*sin(theta_k+x_k(3)-pi/2);0,1,r_k*cos(theta_k+x_k(3)-pi/2)];
 %   g_zk=[cos(theta_k+x_k(3)-pi/2),-r_k*sin(theta_k+x_k(3)-pi/2);sin(theta_k+x_k(3)-pi/2),r_k*cos(theta_k+x_k(3)-pi/2)];


%else
mnew=[r_k*cos(theta_k+x_k(3))+x_k(1);r_k*sin(theta_k+x_k(3))+x_k(2)];
g_xk=[1,0,-r_k*sin(theta_k+x_k(3));0,1,r_k*cos(theta_k+x_k(3))];
g_zk=[cos(theta_k+x_k(3)),-r_k*sin(theta_k+x_k(3));sin(theta_k+x_k(3)),r_k*cos(theta_k+x_k(3))];
%end
   Sigma=[predSigma,predSigma*g_xk';
     g_xk*predSigma, g_xk*predSigma*g_xk'+g_zk*Param.R*g_zk'];
end