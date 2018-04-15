function Li = da_nn(z, Q_t_i)
% perform nearest-neighbor data association

global Param;
global State;

Li = [];

% get robot pose
pos_r = State.Ekf.mu(1:2);
theta = State.Ekf.mu(3);

threshold = chi2inv(0.8,2);%chi2cdf(3,2);
%threshold = chi2inv(0.95,2);
for ind = 1:1:size(z,2)
   r_i = z(1, ind);
   phi_i = z(2, ind);
   pi = [];
   for k = 1:1:State.Ekf.nL
       % get the kth landmark pose(here k is the same order of state vector)
       pos_landmark = State.Ekf.mu((3+2*k-1):(3+2*k));
       delta = pos_landmark - pos_r;
       q = delta' * delta;
       % compute the estimated observation for the kth landmark
       z_hat = [sqrt(q);atan2(delta(2), delta(1))-theta];
       z_hat(2) = minimizedAngle(z_hat(2));
       % compute the projection matrix F_x_k
       F_x_j = zeros(5,length(State.Ekf.mu));
       F_x_j(1:3,1:3) = eye(3);
       F_x_j(4:5,(3+2*k-1):(3+2*k)) = eye(2);
        % compute H_t_i
       H_t_i = 1 / q * [-sqrt(q)*delta(1) -sqrt(q)*delta(2)  0  sqrt(q)*delta(1) sqrt(q)*delta(2);
                         delta(2)         -delta(1)         -q -delta(2)         delta(1)        ] * F_x_j;
       % PSI
       PSI = H_t_i * State.Ekf.Sigma * H_t_i' + Q_t_i;
       delta_z = z(1:2,ind) - z_hat;
       delta_z(2) = minimizedAngle(delta_z(2));
       pi_k = sqrt(delta_z'*pinv(PSI)*delta_z);
       pi = [pi, pi_k];
   end
   [val,index] = min(pi);
   if val <= threshold
       % The smallest distance is smaller than the threshold.
       % index is the index of the landmark
       Li = [Li, index];
   else
       % The smallest distance is larger than the threshold.
       % this is a new landmark
       Li = [Li, 0];
   end
end


end
