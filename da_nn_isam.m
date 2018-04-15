function [Li,measurements] = da_nn_isam(z, Q_t_i)
% perform nearest-neighbor data association

global Param;
global State;

Li = [];
measurements = cell(1, size(z,2));
% get robot pose
pos_r = State.iSAM.rlin(end-2:end-1);
current_pose = pos_r';
theta = State.iSAM.rlin(end);

threshold = chi2inv(0.8,2);%chi2cdf(3,2);

%filter out landmarks by sensor radius
if ~isfield(Param,'sensor_radius')
    Param.sensor_radius=1e6;
end

landmark_list=[];
if State.iSAM.nL>0
    landmark_poss=[State.iSAM.llin(1:2:State.iSAM.nL*2);State.iSAM.llin(2:2:State.iSAM.nL*2)];
    L=sqrt((landmark_poss(1,:)-current_pose(1)).^2+(landmark_poss(2,:)-current_pose(2)).^2)<Param.sensor_radius;
    landmark_list=1:State.iSAM.nL;
    landmark_list= landmark_list(L);
    covariances_extract_naive(landmark_list);
%         covariances_extract(landmark_list);

end

for ind = 1:1:size(z,2)
   r_i = z(1, ind);
   phi_i = z(2, ind);
   pi = [];
   for k = 1:1:length(landmark_list)
       % get the kth landmark pose(here k is the same order of state vector)
       pos_landmark = State.iSAM.llin((2*landmark_list(k)-1):2*landmark_list(k));
       delta = pos_landmark - pos_r;
       q = delta * delta';
       % compute the estimated observation for the kth landmark
       z_hat = [sqrt(q);atan2(delta(2), delta(1))-theta];
       z_hat(2) = minimizedAngle(z_hat(2));
        % compute H_t_i
       H_t_i = 1 / q * [-sqrt(q)*delta(1) -sqrt(q)*delta(2)  0  sqrt(q)*delta(1) sqrt(q)*delta(2);
                         delta(2)         -delta(1)         -q -delta(2)         delta(1)        ];
       % PSI
       predSigma = State.iSAM.sigma([1:3,3+2*k-1,3+2*k],[1:3,3+2*k-1,3+2*k]);
       PSI = H_t_i * predSigma * H_t_i' + Q_t_i;
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
       measurements{ind} = [];
   else
       % The smallest distance is larger than the threshold.
       % this is a new landmark
       Li = [Li, 0];
       measurements{ind} = [State.iSAM.nR;z(1:3,ind);State.iSAM.rlin(State.iSAM.nR*3-2:State.iSAM.nR*3)'];
   end
end


end
