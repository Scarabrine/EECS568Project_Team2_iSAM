function [Li] = da_known_pre(z,nodes)
% Here it is a hack version of known data association
global Param;
global State;
Li = [];
if nodes(2) > State.iSAM.max_node
    %% This observation is a new node
    Li = [0];
else
    %% This observation is not a new node
    index = find(State.iSAM.lNode == nodes(2));
    if length(index) ==1
        Li = [index];
    else
        error('Observed landmark node should only be in the lNode once and only once.')
    end
end
end

