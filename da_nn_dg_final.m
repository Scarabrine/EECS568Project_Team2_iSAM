function [Li,measurements]=da_nn_dg_final(z)
% perform nearest-neighbor data association
%need to
%update for isam to get covariance matrix we need
%modify to spit out cell called measurements of size 1xNobs. Include of pose number, observation, for each one
%registered with landmark for intiialization. Leave cell empty if obs is not being used or landmark already initialized
%Li=1xNobs array with landmark id. 0 if new landmark, -1 if skip

global Param;
global State;

%get size of observations and preallocate index list
m=size(z,2);
measurements=cell(1,m);
Li=-1*ones(1,m);
threshold=chi2inv(0.9 ,2); %0.95
nHits_threshold=10;%number of hits required before landmark is initialized
%number of hits required before landmark is initialized
%not sure if this is right or if we need to add anything to the
%linearization point

if strcmp(Param.method,'iSAM')
   %landmark positions and covariance
   landmark_positions=[State.iSAM.llin(1:2:State.iSAM.nL*2);State.iSAM.llin(2:2:State.iSAM.nL*2)];
   covariances_extract_naive(1:State.iSAM.nL);
   landmark_covariances=cell(1,State.iSAM.nL);
   for k=1:State.iSAM.nL
       landmark_covariances{k}=State.iSAM.sigma([1:3,3+2*k-1,3+2*k],[1:3,3+2*k-1,3+2*k]);
   end
   %Use sensor model to get zhat
   current_pose(1,1) = State.iSAM.rlin(end-2);
   current_pose(2,1) = State.iSAM.rlin(end-1);
   current_pose(3,1)= State.iSAM.rlin(end);
   tstep=State.iSAM.nR;
elseif strcmp(Param.method,'Ekf')
   %Index state vector to get x,y of the landmark location
   landmark_positions = [State.Ekf.mu(4:2:3+2*State.Ekf.nL)';State.Ekf.mu(5:2:3+2*State.Ekf.nL)'];
   landmark_covariances=cell(1,State.Ekf.nL);
   for k=1:State.Ekf.nL
       landmark_covariances{k}=State.Ekf.Sigma([1:3,3+2*k-1,3+2*k],[1:3,3+2*k-1,3+2*k]);
   end
   %Use sensor model to get zhat
   current_pose(1,1) = State.Ekf.mu(1);
   current_pose(2,1) = State.Ekf.mu(2);
   current_pose(3,1) = State.Ekf.mu(3); 
   tstep=State.Ekf.t;
end

%reduce threshold if we dont have any landmarks yet
if ~exist('landmark_positions','var')
    nHits_threshold=round(nHits_threshold/2);
end

%initialize structure to hold proposed landmarks
if ~isfield(State,'proposed_nL')
    State.proposed_nL=0;
    State.mu_proposed=[];
    State.nHits=[];
    State.proposed_average_measurement=cell([0 0]);
    State.proposed_individual_measurements=cell([0 0]);
    State.proposed_landmark_covariances=cell([0 0]);
end
% construct covariances for proposed landmarks
      % for updating proposed entry
      if strcmp(Param.method,'iSAM')
    stateSigma=State.iSAM.sigma(1:3,1:3);
elseif strcmp(Param.method,'Ekf')
       stateSigma=State.Ekf.Sigma(1:3,1:3);
        end

%nearest neighbor algorithm for list of proposed and current landmarks
for i=1:m
  [Lindex,D2]=nn_compatibility(z(:,i),current_pose,landmark_positions,landmark_covariances);
  if D2<threshold
      Li(i)=Lindex;
      measurements{i}=[];
  else
  [Lindex,D2]=nn_compatibility(z(:,i),current_pose,State.mu_proposed,State.proposed_landmark_covariances);
  [mnew,Signew]=get_pred_landmark_position(z(:,i),current_pose,stateSigma);
  if D2<threshold
      State.mu_proposed(:,Lindex)=(State.mu_proposed(:,Lindex)*State.nHits(Lindex)+mnew)/(State.nHits(Lindex)+1);
      State.proposed_landmark_covariances{Lindex}=(State.proposed_landmark_covariances{Lindex}*State.nHits(Lindex)+Signew)/(State.nHits(Lindex)+1);
      State.nHits(Lindex)= State.nHits(Lindex)+1;
      State.proposed_individual_measurements{Lindex}=[State.proposed_individual_measurements{Lindex},[tstep;z(1:3,i);current_pose]];
      State.proposed_average_measurement{Lindex}=[tstep;get_measurement_average(State.mu_proposed(:,Lindex),current_pose)];
      if State.nHits(Lindex)>nHits_threshold
          %tell to initialize new landmark
          Li(i)=0;
          if strcmp(Param.method,'iSAM')
              %measurements{i} = [State.iSAM.nR;z(1:3,i);State.iSAM.rlin(State.iSAM.nR*3-2:State.iSAM.nR*3)'];
              measurements{i}=State.proposed_average_measurement{Lindex};
              %measurements{i}=State.proposed_individual_measurements{Lindex};
              
          else
         measurements{i}=State.proposed_average_measurement{Lindex};
         %measurements{i} = [tstep;z(1:3,i);current_pose];
          end
          %delete landmark from proposed list
          State.proposed_nL=State.proposed_nL-1;
          State.nHits(Lindex)=[];
          State.mu_proposed(:,Lindex)=[];
          State.proposed_average_measurement(Lindex)=[];
          State.proposed_individual_measurements(Lindex)=[];
          State.proposed_landmark_covariances(Lindex)=[];
      end
  else
      %add to proposal list
      measurements{i}=[];
      State.proposed_nL=State.proposed_nL+1;
      State.nHits(State.proposed_nL)=1;
      State.mu_proposed(:,State.proposed_nL)=mnew;
      State.proposed_average_measurement{State.proposed_nL}=z(1:2,i);
      State.proposed_individual_measurements{State.proposed_nL}=[tstep;z(1:3,i);current_pose];
      State.proposed_landmark_covariances{State.proposed_nL}=Signew;
  end
end

end

