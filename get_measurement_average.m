function [z_hat]=get_measurement_average(mu,current_pose)
 delta = mu - current_pose(1:2);
       q = (delta')* delta;
       z_hat = [sqrt(q);atan2(delta(2), delta(1))-current_pose(3)];
        z_hat(2)=minimizedAngle(z_hat(2));
end