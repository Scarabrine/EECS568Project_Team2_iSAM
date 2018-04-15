function  recover (i,l,inv_R)

global complete_sigma
global State

sum = 0;
[~,l2]=size(complete_sigma);
for j=find(State.iSAM.R(i,:))
    if (j~=i)
        if (j>l)
          if (isnan(complete_sigma(l,j)))
              recover(l, j, inv_R);
          end
          sum = sum + State.iSAM.R(i,j) * complete_sigma(l,j);
        else
          if (isnan(complete_sigma(j,l)))
              recover(j, l, inv_R);
          end
          sum = sum + State.iSAM.R(i,j) * complete_sigma(j,l);
        end
    end
end

if (i == l)                                 % diagonal entries
    complete_sigma(i,l) =  (inv_R(l) * (inv_R(l) - sum));
    if(isnan(complete_sigma(i,l)))
        error('error in recover function. Wrongly set sigma value as nan');
    end
else                                        % off-diagonal entries
    complete_sigma(i,l) = (- sum * inv_R(i));
    if(isnan(complete_sigma(i,l)))
        error('error in recover function. Wrongly set sigma value as nan');
    end
    complete_sigma(l,i) = complete_sigma(i,l);
end