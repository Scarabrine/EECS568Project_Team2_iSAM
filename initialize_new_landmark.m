function initialize_new_landmark(z,nodes)

global Param;
global State;

State.iSAM.max_node = State.iSAM.max_node+1;
State.iSAM.lNode = [State.iSAM.lNode, nodes(2)];

c = cos(State.Ekf.mu(3)); s = sin(State.Ekf.mu(3)); 


m_x = z(1)*c - z(2)*s + State.Ekf.mu(1);
m_y = z(1)*s + z(2)*c + State.Ekf.mu(2);

dx = m_x - State.Ekf.mu(1); 
dy = m_y - State.Ekf.mu(2);

State.Ekf.mu = [State.Ekf.mu; m_x; m_y];
d= sqrt(z(1)^2+z(2)^2);


J1 = [ c , s;
      -s , c];
  
J2 = [-c, -s, -s*dx + c*dy;
      s, -c, -c*dx - s*dy];
    
    
State.Ekf.Sigma = [             State.Ekf.Sigma        ,         zeros(length(State.Ekf.Sigma),2)      ;
                    zeros(2,length(State.Ekf.Sigma))   ,  J2*State.Ekf.Sigma(1:3,1:3)*J2' + J1*Param.R*J1'  ];

State.Ekf.Sigma(end-1:end,1:3) = J2*State.Ekf.Sigma(1:3,1:3);
State.Ekf.Sigma(1:3,end-1:end) = (State.Ekf.Sigma(end-1:end,1:3))';

