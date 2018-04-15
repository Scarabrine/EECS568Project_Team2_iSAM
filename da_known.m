function [Li,measurements] = da_known(z, sL)
% data association with known correspondences
% z: observation, containts several landmark observations
% sL: a list of observed landmark signatures. sL(i)=j means the ith
% landmark is named with signature j.
% Li: a vector with the size of observations. Li(i) = j means that the ith
% observation is landmark j. If j=0, then it is a new landmark.


global Param;
global State;

Li = [];
measurements = cell(1, size(z,2));
% For each observation, check the correspondence index of this landmark in
% the state vector
% 0. Get the observed landmark ids(signatures)
observed_landmark = sL;
for ind = 1:1:size(z,2)
    % 1. Get the current landmark signature for current observation
    signature = z(3,ind);
    % 2. Check if the signature is in the observed_landmark list or not. If
    % not set to be 0. otherwise return the index
    current_index = find(observed_landmark == signature);
    if(~isempty(current_index))
        % 2.1 this landmark has been observed before. 
        Li = [Li, current_index];
        measurements{ind} = [];

    else
        % 2.2 this landmark is new
        Li = [Li, 0];
        measurements{ind} = [State.iSAM.nR;z(1:3,ind);State.iSAM.rlin(State.iSAM.nR*3-2:State.iSAM.nR*3)'];

    end
end

end