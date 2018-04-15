function covariances_extract_naive (landmark_nos)

% input = array of landmark ids in order
% output = square covariance matrix of robot pose followed by landmarks in
%          same order as the input

global State

% indices in R of current robot x,y,theta and required landmarks
index_robot = [State.iSAM.rM(3*State.iSAM.nR-2), State.iSAM.rM(3*State.iSAM.nR-1), State.iSAM.rM(3*State.iSAM.nR)];
index_landmark=[];
for i=1:size(landmark_nos,2)
    index_landmark = [index_landmark, State.iSAM.lM(2*landmark_nos(i)-1),State.iSAM.lM(2*landmark_nos(i))];
end
indices =[index_robot,index_landmark];
comp_sig = inv(State.iSAM.R'*State.iSAM.R);
State.iSAM.sigma = comp_sig(indices,indices);

end
  