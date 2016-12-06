function r1 = rVsFlipPercV2(obj,nSampleIntervals)

%*************************************************************************
% This function is used to calculate the R factor over the flipping
% percentage. Note that only lowest charges under the percentage are flipped.
% Created by Yun Zhao on May 20th, 2015
%
% Comment: It can be easily editted to calculate R factor over absolute 
% charge densities.
% Compared with first Version: Here use chargeFlipPercentV2 to flip.
%*************************************************************************

flipPerc = linspace(0,100,nSampleIntervals)/100;
r1 = zeros(nSampleIntervals,2);
mask = ones(size(obj)); % Because mask is a required input for function getRfactor
Fobj = abs(fftn(obj));

for i = 1:nSampleIntervals
    objFlip = chargeFlipPercentV2(obj,flipPerc(i));
    FobjFlip = abs(fftn(objFlip));
    r1(i,:) = [i/nSampleIntervals getRFactor(Fobj,FobjFlip,mask)];
end

figure
plot(r1(:,1),r1(:,2));
xlabel('flipping fraction of charge values in real space','fontsize',14);
ylabel('R factor','fontsize',14);
title('charge flipping effects on R factor','fontsize',14);

    
