function initialize_new_landmark(z)

global Param;
global State;

%get predicted pose
x_k=State.Ekf.mu(1:3);
%predicted measurement 
r_k=z(1);
theta_k=z(2);

%estimate new landmark position

mnew=[r_k*cos(theta_k+x_k(3))+x_k(1);r_k*sin(theta_k+x_k(3))+x_k(2)];

%augment state
State.Ekf.mu=[State.Ekf.mu;mnew];

%augment covariance matrix
g_xk=[1,0,-r_k*sin(theta_k+x_k(3));0,1,r_k*cos(theta_k+x_k(3))];

g_zk=[cos(theta_k+x_k(3)),-r_k*sin(theta_k+x_k(3));sin(theta_k+x_k(3)),r_k*cos(theta_k+x_k(3))];

if State.Ekf.nL==0
    State.Ekf.Sigma=[State.Ekf.Sigma,State.Ekf.Sigma*g_xk';
      g_xk*State.Ekf.Sigma, g_xk*State.Ekf.Sigma*g_xk'+g_zk*Param.R*g_zk'];
else
    R=Param.R;
    Pxx=State.Ekf.Sigma(State.Ekf.iR,State.Ekf.iR);
    Pmm=State.Ekf.Sigma(State.Ekf.iM,State.Ekf.iM);
    Pxm=State.Ekf.Sigma(State.Ekf.iR,State.Ekf.iM);
    State.Ekf.Sigma=[Pxx,Pxm,Pxx'*g_xk';...
                    Pxm',Pmm,Pxm'*g_xk';...
                    g_xk*Pxx,g_xk*Pxm,g_xk*Pxx*g_xk'+g_zk*R*g_zk'];
end

%record indices
State.Ekf.nL=State.Ekf.nL+1;
State.Ekf.sL=[State.Ekf.sL,State.Ekf.nL];
State.Ekf.iL{State.Ekf.nL}=[length(State.Ekf.mu)-1,length(State.Ekf.mu)];
State.Ekf.iM=State.Ekf.iR(end)+1:length(State.Ekf.mu);
end