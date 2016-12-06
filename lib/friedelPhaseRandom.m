function PhaseFridRand = friedelPhaseRandom(phaseSize,flag)
%********************************************************************
%Creat an random phase array to an array which satisfy Friedel's law
%*********************************************************************

phase = rand(phaseSize).*2*pi;
dimNum = length(phaseSize);
switch dimNum
    case 1
        objFlip = flip(phase,1);
    case 2
        objFlip = flip(flip(phase,1),2);
    case 3
        objFlip = flip(flip(flip(phase,1),2),3);
    otherwise
        fprintf('Invalid dimension\n');
end


PhaseFridRand1 = phase - objFlip;
% The center of phi(000) is in the middle of array "PhaseFridRand" which is
% 0. And I always assume the dimension sizes are odd


switch flag
    case 'centered'
        PhaseFridRand = PhaseFridRand1;
    case 'cornered'
        PhaseFridRand = ifftshift(PhaseFridRand1);
    otherwise
        PhaseFridRand = ifftshift(PhaseFridRand1);
end

