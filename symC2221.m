function objSymAve = symC2221(triObj)

%% Note 
% In function symC2221_v2.m, it assumes that input obj is one unit cell.
% In function symC2221.m, it assumes that input is a triple cell. zeros are
% padded above and below the unit cell.
%% This function impose symmetry constraint from real space.
% Symmetry operation only for space group C 2 2 21.
% Example from 3rdx.pdb
% REMARK 290      SYMOP   SYMMETRY                                                
% REMARK 290     NNNMMM   OPERATOR                                                
% REMARK 290       1555   X,Y,Z                                                   
% REMARK 290       2555   -X,-Y,Z+1/2                                             
% REMARK 290       3555   -X,Y,-Z+1/2                                             
% REMARK 290       4555   X,-Y,-Z                                                 
% REMARK 290       5555   X+1/2,Y+1/2,Z                                           
% REMARK 290       6555   -X+1/2,-Y+1/2,Z+1/2                                     
% REMARK 290       7555   -X+1/2,Y+1/2,-Z+1/2                                     
% REMARK 290       8555   X+1/2,-Y+1/2,-Z                                         

%% Matrix after symmetry operation and averaging.

sizeC = size(triObj,3);
obj = triObj(:,:,int32(sizeC/3+1):int32(sizeC*2/3)); 


objSize = size(obj);


xShift = (objSize(1)-1)/2;
yShift = (objSize(2)-1)/2;
zShift = (objSize(3)-1)/2;

obj2 = circshift(flip(flip(obj,1),2),zShift,3);
obj3 = circshift(flip(flip(obj,1),3),zShift,3);
obj4 = flip(flip(obj,2),3);
obj5 = circshift(circshift(obj,xShift,1),yShift,2);
obj6 = circshift(circshift(circshift(flip(flip(obj,1),2),xShift,1),yShift,2),zShift,3);
obj7 = circshift(circshift(circshift(flip(flip(obj,1),3),xShift,1),yShift,2),zShift,3);
obj8 = circshift(circshift(flip(flip(obj,2),3),xShift,1),yShift,2);


objSymAve1 = (obj+obj2+obj3+obj4+obj5+obj6+obj7+obj8)/8;


%% Padd to a triple cell.
objSymAve = zeros(size(triObj));
objSymAve(:,:,int32(sizeC/3+1):int32(sizeC*2/3)) = objSymAve1;




%% operation matrix. (My initial idea. not good. I keep it here just for memory)
% Yun's note: I just realized that we don't need to do this matrix
% rotation. The new coordinate after symmetry operation is already shown in
% REMARK. We can actually use flip and circshift function in matlab to
% implement these symmetry operations. After that, we avearge them. It
% should be very fast.

% sym1 = [1 0 0;0 1 0;0 0 1];
% sym2 = [-1 0 0;0 -1 0;0 0 1];
% sym3 = [-1 0 0;0 1 0;0 0 -1];
% sym4 = [1 0 0;0 -1 0;0 0 -1];
% sym5 = [1 0 0;0 1 0;0 0 1];
% sym6 = [-1 0 0;0 -1 0;0 0 1];
% sym7 = [-1 0 0;0 1 0;0 0 -1];
% sym8 = [1 0 0;0 -1 0;0 0 -1];
% 
% r1 = [0 0 0];
% r2 = [0 0 0.5];
% r3 = [0 0 0.5];
% r4 = [0 0 0];
% r5 = [0.5 0.5 0];
% r6 = [0.5 0.5 0.5];
% r7 = [0.5 0.5 0.5];
% r8 = [0.5 0.5 0];


