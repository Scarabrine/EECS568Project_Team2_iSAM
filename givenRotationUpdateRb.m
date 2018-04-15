function [] = givenRotationUpdateRb(k_start,k_end,i_start, i_end)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global State
[m,n] = size(State.iSAM.R);
%% Parameter
eps = 1e-10;
for k = k_start:1:k_end
    for i = max(k+1,i_start):i_end
        if State.iSAM.R(i,k) == 0
            continue;
        end
        
        % Given rotation matrix
        alpha = State.iSAM.R(k,k);
        beta = State.iSAM.R(i,k);
        
        if beta == 0
            c = 1;
            s = 0;
        elseif abs(beta) > abs(alpha)
            c = -alpha/beta/sqrt(1+(alpha/beta)^2);
            s = 1/sqrt(1+(alpha/beta)^2);
        else
            c = 1/sqrt(1+(beta/alpha)^2);
            s = -beta/alpha/sqrt(1+(beta/alpha)^2);
        end
        G = speye(m);
        G(k,k) = c;
        G(k,i) = s;
        G(i,k) = -s;
        G(i,i) = c;
         
        State.iSAM.R = G' * sparse(State.iSAM.R);%% update R
        State.iSAM.b = G' * sparse(State.iSAM.b);%% update b
        if abs(State.iSAM.R(i,k)) < abs(eps)
            State.iSAM.R(i,k) = 0;
        end
        
    end
end

end

