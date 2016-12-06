function moduleFridRand = friedelModuleRandom(moduleAverage,moduleSize)
%*********************************************************
%Convert a modulus array to an array which satisfy Friedel's law
%*********************************************************

module = rand(moduleSize).*moduleAverage;
objFlip = flip(flip(flip(module,1),2),3);
moduleFridRand = (module + objFlip)/2;
% The diagonal elements of obj2 should be zero if that dimension size is odd.
