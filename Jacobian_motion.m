function [J_state,J_control] = Jacobian_motion (u)

global State;
state = State.Ekf.mu(1:3);

J_state = eye(3);
J_state (1:2,3) = [-u(1)*sin(state(3))-u(2)*cos(state(3)); u(1)*cos(state(3))-u(2)*sin(state(3))];

J_control = [cos(state(3))  ,   -sin(state(3))   ,   0;
             sin(state(3)) ,   cos(state(3))   ,   0;
                    0       ,         0         ,    1    ];

end
