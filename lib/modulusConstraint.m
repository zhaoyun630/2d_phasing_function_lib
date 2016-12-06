function modulusNew = modulusConstraint(FObserved,FEstimated,mask)

%**********************************************************************
%       Used to enforce modulus constraint with a given mask. 
%       Use this function when only part of structure factor amplitudes need
% to be replaced with measured ones.
%**********************************************************************
if mask == 0
    % 0 means no mask is provided.
        modulusNew = FObserved;
else
        modulusNew = FObserved.*mask + FEstimated.*(1-mask);
end
            
