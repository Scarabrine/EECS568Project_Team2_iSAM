%Returns true if hypothesis and measurements are jointly compatible
%Measurements is 2 x N, N is number of measurements.
%Hypothesis is 1 x N. Is indices for each guess (hypothesis) of each
%measurement
function result = joint_compatibility(measurements,hypothesis)
    global Param;
    global State;
    global JCBB;
    nummeasurements = size(measurements,2);
    
    %Maithili's function. Get covariance
    %covariances_extract(hypothesis)
    %cov = State.iSAM.sigma;
    %UPDATE: Don't call Maithili's function. Too slow.
    stateindices = [1 2 3];
    landmarkindices = getLandmarkIndices(hypothesis);
    desiredindices = [stateindices landmarkindices];
    cov = JCBB.cov(desiredindices,desiredindices);
    
    R = [];
    H = [];
    differences = [];
    %Create stacked H matrix. 1 H matrix for each landmark in hypothesis
    for i=1:nummeasurements
       obs = measurements(1:2,i);
       landmarkid = hypothesis(:,i); 
       
       if strcmp(Param.method,'iSAM')
           %Index state vector to get x,y of the landmark location
           lm_locationx = State.iSAM.llin(2 * landmarkid -1);
           lm_locationy = State.iSAM.llin(2 * landmarkid);       
           lm_location = [lm_locationx;lm_locationy];

           %Use sensor model to get zhat
           robotx = State.iSAM.rlin(end-2);
           roboty = State.iSAM.rlin(end-1);
           robottheta= State.iSAM.rlin(end);
       elseif strcmp(Param.method,'EKF')
           %Index state vector to get x,y of the landmark location
           lm_locationx = State.Ekf.mu(3 + 2 * landmarkid -1);
           lm_locationy = State.Ekf.mu(3 + 2 * landmarkid);       
           lm_location = [lm_locationx;lm_locationy];

           %Use sensor model to get zhat
           robotx = State.Ekf.mu(1);
           roboty = State.Ekf.mu(2);
           robottheta = State.Ekf.mu(3);           
       end

       delta = [lm_location(1)-robotx;
                lm_location(2)-roboty];
       q = delta' * delta;

%       if strcmp(Param.choice,'vp')
% Victoria park has a different sensor model???
%           temp = minimizedAngle(atan2(delta(2),delta(1)) - robottheta + pi/2);
%       else
           temp = minimizedAngle(atan2(delta(2),delta(1)) - robottheta);
%       end


       %Sensor model to get zhat
       zhat = [sqrt(q)
                temp];
       %Compute F matrix
       %TODO: Fix j value. Since covariance was truncated?
       j = i;
       %Truncate number of landmarks
       numlandmarks = nummeasurements;
       Fxj = [eye(3)       zeros(3, 2 * j - 2)   zeros(3,2)  zeros(3,2 * numlandmarks - 2 * j);
              zeros(2,3)   zeros(2, 2 * j - 2)   eye(2)      zeros(2,2 * numlandmarks - 2 * j)];

       %Get Jacobian
       sq = sqrt(q);   
       Hit = (1/q) * [-sq*delta(1) -sq*delta(2)  0  sq*delta(1) sq*delta(2);
                       delta(2)    -delta(1) -q    -delta(2)    delta(1)];
       Hit = Hit * Fxj;
       
      temp = obs - zhat;
      temp(2) = minimizedAngle(temp(2));
        
       differences = vertcat(differences,temp);
       H = vertcat(H,Hit);
       R = blkdiag(R,JCBB.R);
   end
        
   %Get covariance for mahalanobis distance
   Psi = H * cov * H' + R;
   
   %calculate mahalanobis distance
   chi2distance = (differences)' * Psi^-1 * (differences);
   
   dof = nummeasurements * 2;
   
   result = false;
   threshold = chi2inv(1-JCBB.alpha,dof);
   if chi2distance < threshold
       result = true;
   end
      
end

function landmarkindices = getLandmarkIndices(landmarkindices)
    lx = 3 + (landmarkindices * 2 - 1);
    ly = 3 + (landmarkindices * 2);
    
    %Interleave lx and ly
    A = [lx;ly];   % concatenate them vertically
    landmarkindices = A(:); %Flatten the result
    landmarkindices = reshape(landmarkindices,1,numel(landmarkindices));
end