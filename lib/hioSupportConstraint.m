function [objNew,rms] = hioSupportConstraint(objK,objK1,support,feedbackParameter,errorflag)
%**************************************************************************
%   Enforce real space constraint from support in HIO algorithm.
%   Note:objK is from last iteration.
%           objK1 is the inverse Fourier transform of structure factors in 
%           current iteration.
%**************************************************************************
objNew = abs(objK1.*support + objK*(1-support)-feedbackParameter.*objK1(1-support));

if errorflag=0
	rms =[];
else
	densityOutsideSupport = abs(objK1).*(1-support);
	densityTotal = abs(objK1);
	rms = sqrt(sum(sum(sum(densityOutsideSupport.^2)))/(sum(sum(sum(densityTotal.^2)))));
end





