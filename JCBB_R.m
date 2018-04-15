%-------------------------------------------------------
function JCBB_R (observations, H, i)
% 
%-------------------------------------------------------
global JCBB;
global State;
global Param;

%i represents the measurement number
if i > JCBB.m % leaf node?
    if pairings(H) > pairings(JCBB.Best) % did better?
        JCBB.Best = H;
    end
else
    if strcmp(Param.method,'iSAM')
        numlandmarks = State.iSAM.nL;
    elseif strcmp(Param.method,'EKF')
        numlandmarks = State.Ekf.nL;
    end
    
    %Loop through all landmarks
    for j=1:numlandmarks
        %Check individual compatibility
        if individual_compatibility(JCBB.z(:,i),[j])
            %Check joint compatibility. We have z...and H
            H_aug = [H j];
            indices = find(H_aug);
            if joint_compatibility(JCBB.z(:,indices),H_aug(:,indices))
                JCBB_R(observations, [H j], i + 1); %pairing (Ei, Fj) accepted 
            end
        end
    end
    
%     individually_compatible = find(compatibility.ic(i,:));
%     for j = individually_compatible
%         if jointly_compatible(prediction, observations, [H j])
%             JCBB_R(observations, compatibility, [H j], i + 1); %pairing (Ei, Fj) accepted 
%         end
%     end
    
%    if pairings(H) + length(compatibility.candidates.observations) - i >= pairings(Best) % can do better?
     if pairings(H) + JCBB.m - i >= pairings(JCBB.Best) % can do better?
%    if pairings(H) + pairings(compatibility.AL(i+1:end)) >= pairings(Best) % can do better?
        JCBB_R(observations, [H 0], i + 1); % star node: Ei not paired
     end
end

end

%-------------------------------------------------------
% 
%-------------------------------------------------------
%Gets the number of nonzero hypotheses
function p = pairings(H)
p = length(find(H));
end
