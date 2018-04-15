%z is a 2 x number of observations matrix
%R is sensor noise covariance matrix. should be 2x2.
%Returns association, a 1 x (number of observations) matrix. Each index
%represents the best landmark assignment.
function [association,measurements] = da_jcbb(z, R)
% perform joint-compatability branch and bound data association

global JCBB;
global State;
global Param;

%Tunable Param. Depends on what Mahalanobis distance you want
JCBB.alpha = .001;

JCBB.R = R; %Sensor noise. Range and bearing
JCBB.m = size(z,2); %Number of measurements
JCBB.z = z; %Measurements

JCBB.Best = zeros(1, size(z,2));

%Extract covariance matrix once. Helps so that we don't need to keep
%recomputing it
%Get all of the covariance matrix. covariances_extract will give robot pose
%covariance

if strcmp(Param.method,'iSAM')
    matdesired = 1:State.iSAM.nL;
    covariances_extract(matdesired);

    %Use this and continuously rearrange it to get desired covariance elements
    JCBB.cov = State.iSAM.sigma;
elseif strcmp(Param.method,'EKF')
    JCBB.cov = State.Ekf.Sigma;
end

JCBB_R (z, [], 1);

association = JCBB.Best;
measurements = cell(1,JCBB.m);
indices = 1:size(z,2);
for ind = indices
    if strcmp(Param.method,'iSAM')
        measurements{ind} = [State.iSAM.nR;z(1:3,ind);State.iSAM.rlin(State.iSAM.nR*3-2:State.iSAM.nR*3)'];
    elseif strcmp(Param.method,'EKF')
        %Added association(ind)
        measurements{ind} = [State.Ekf.t;z(1:2,ind); association(ind); State.Ekf.mu(1:3,1)];
    end
end

%Change to -1 where we want new landmarks
%association(association == 0) = -1;
end

